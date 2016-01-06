#include <iostream>
#include "Img.h"
#include <netcdfcpp.h>
// to compile with netcdfcpp.h, add -lnetcdf_c++
#include "Spore.h"
//g++ -std=c++11 main.cpp Img.h Img.cpp Spore.h Spore.cpp -lgdal -lnetcdf_c++


using namespace std;

#define UMCA_RAST_PATH "./layers/UMCA_den_100m.img"
#define OAKS_RAST_PATH "./layers/OAKS_den_100m.img"
#define LVTREE_RAST_PATH "./layers/TPH_den_100m.img"
#define I_OAKS_RAST_PATH "./layers/init_2000_cnt.img"
#define BACKGROUND_IMAGE_PATH "./layers/ortho_5m_color.tif"
#define START_TIME 2000
#define END_TIME 2010
#define WEATHER_COEFF_PATH "./layers/weather/weatherCoeff_2000_2014.nc"
#define DIM 1

// Initialize infected trees for each species(!!Needed unless empirical info is available)
static Img initialize(Img& img1,Img& img2){
	int re_width=0;
	int re_height=0;
	int **re_data=NULL;
	if(img1.getWidth() != img2.getWidth() || img2.getHeight() != img2.getHeight()){
		cerr << "The height or width of one image do not match with that of the other one!" << endl;
		return Img();
	}else{
		re_width = img1.getWidth();
		re_height = img1.getHeight();
		re_data = (int **)CPLMalloc(sizeof(int *) * re_height);
		int *stream = (int *)CPLMalloc(sizeof(int)*re_width*re_height);

		for(int i=0;i<re_height;i++){
			re_data[i] = &stream[i*re_width];
		}

		for(int i=0;i<re_height;i++){
			for(int j=0;j<re_width;j++){
				if(img2.data[i][j]>0){
					if(img1.data[i][j]>img2.data[i][j])
						re_data[i][j] = img1.data[i][j] < (img2.data[i][j]*2)?img1.data[i][j]:(img2.data[i][j]*2);
					else
						re_data[i][j] = img1.data[i][j];
				}else{
					re_data[i][j]=0;
				}
			}
		}

		return Img(re_width,re_height,img1.getWEResolution(), img1.getNSResolution(),re_data);
	}
}


int main()
{
	// set the start point of the program
	clock_t begin = clock();

	// read the suspectible UMCA raster image
	Img umca_rast(UMCA_RAST_PATH);

	// read the SOD-affected oaks raster image
	Img oaks_rast(OAKS_RAST_PATH);

	// read the living trees raster image
	Img lvtree_rast(LVTREE_RAST_PATH);

	// read the initial infected oaks image
	Img I_oaks_rast(I_OAKS_RAST_PATH);

	// create the initial suspectible oaks image
	Img S_oaks_rast = oaks_rast - I_oaks_rast;

	// create the initial infected umca image
	Img I_umca_rast = initialize(umca_rast,I_oaks_rast);

	// create the initial suspectible umca image
	Img S_umca_rast = umca_rast - I_umca_rast;

	// SOD-immune trees image
	//Img SOD_rast = umca_rast + oaks_rast;
	//Img IMM_rast = lvtree_rast - SOD_rast;

/*
	for(int i=0;i<S_oaks_rast.getHeight();i++){
		for(int j=0;j<S_oaks_rast.getWidth();j++){
			cout << S_oaks_rast.data[i][j] << " ";
		}
		cout << endl;
	}
*/
	// read the background satellite image
	Img bkr_img(BACKGROUND_IMAGE_PATH);

	// retrieve the width and height of the images
	int width = umca_rast.getWidth();
	int height = umca_rast.getHeight();

	if(START_TIME > END_TIME){
		cerr << "Start date must precede the end date!!!"<< endl;
		exit(EXIT_FAILURE);
	}

	NcFile weatherCoeff(WEATHER_COEFF_PATH, NcFile::ReadOnly);

	if(!weatherCoeff.is_valid()){
		cerr << "Can not open the weather coefficients file(.cn)!"<< endl;
		exit(EXIT_FAILURE);
	}

	NcVar *mcf_nc, *ccf_nc;
	if(!(mcf_nc = weatherCoeff.get_var("Mcoef"))){
		cerr << "Can not read the moisture coefficients from the file!"<< endl;
		exit(EXIT_FAILURE);
	}

	if(!(ccf_nc = weatherCoeff.get_var("Ccoef"))){
		cerr << "Can not read the temperature coefficients from the file!"<< endl;
		exit(EXIT_FAILURE);
	}

	double mcf[height][width];
	double ccf[height][width];
	double **weather = new double *[height];
	for(int i=0;i<height;i++){
		weather[i] = new double[width];
	}

	// Seasonality: Do you want the spread to be limited to certain months?
	bool ss = true;
	bool wind = true;
	// set the direction of the wind {N=0,NE=45,E=90,SE=135,S=180,SW=225,W=270,NW=315,NO  // NO means that there is no wind}
	Direction pwdir=NE;
	// set the spore rate
	double spore_rate = 4.4;
	Rtype rtype = CAUCHY;
	double scale1 = 20.57;
	int kappa = 2;

	// initialize the start Date and end Date object
	Date dd_start(START_TIME,01,01);
	Date dd_end(END_TIME,12,31);

	// main simulation loop(weekly steps)
	for(int i=0;dd_start.compareDate(dd_end);i++,dd_start.increasedByWeek()){
		bool allInfected = true;
		for(int j=0;j<height;j++){
			for(int k=0;k<width;k++){
				if(S_oaks_rast.data[j][k]>0)
					allInfected = false;
			}
		}

		// if all the oaks are infected, then exit
		if(allInfected){
			cerr << "In the " << dd_start.getYear() << "-" << dd_start.getMonth() <<"-" << dd_start.getDay() << " all suspectible oaks are infected!" << endl;
			break;
		}

		// check whether the spore occurs in the month
		if(ss && dd_start.getMonth()>9){
			cout << "----------"  << dd_start.getYear() << "-" << dd_start.getMonth() <<"-" << dd_start.getDay() << "-------------" << endl;
			for(int m=0;m<height;m++){
				for(int j=0;j<width;j++){
					cout << I_oaks_rast.data[m][j] << " ";
				}
				cout << endl;
			}
			continue;
		}
		
		// read the weather information
		if(!mcf_nc->set_cur(i,0,0)){
			cerr << "Can not read the coefficients from the mcf_nc pointer to mcf array " << i << endl;
			exit(EXIT_FAILURE);
		}

		if(!ccf_nc->set_cur(i,0,0)){
			cerr << "Can not read the coefficients from the ccf_nc pointer to ccf array "<< i << endl;
			exit(EXIT_FAILURE);
		}

		if(!mcf_nc->get(&mcf[0][0],1,height,width)){
			cerr << "Can not get the record from mcf_nc "<< i << endl;
			exit(EXIT_FAILURE);
		}

		if(!ccf_nc->get(&ccf[0][0],1,height,width)){
			cerr << "Can not get the record from ccf_nc "<< i << endl;
			exit(EXIT_FAILURE);
		}

		for(int j=0;j<height;j++){
			for(int k=0;k<width;k++){
				weather[j][k] = mcf[j][k] * ccf[j][k];
			}
		}

		// build the Sporulation object
		Sporulation sp1;
		sp1.SporeGen(I_umca_rast,weather,spore_rate);

		sp1.SporeSpreadDisp(S_umca_rast, S_oaks_rast, I_umca_rast, I_oaks_rast, lvtree_rast, 
				rtype, weather, scale1, kappa, pwdir);

		// print out the result
		cout << "----------"  << dd_start.getYear() << "-" << dd_start.getMonth() <<"-" << dd_start.getDay() << "-------------" << endl;
		for(int m=0;m<height;m++){
			for(int j=0;j<width;j++){
				cout << I_oaks_rast.data[m][j] << " ";
			}
			cout << endl;
		}

	}

	// compute the time used when running the model
	clock_t end = clock();
	double elapsed_secs = double(end-begin) / CLOCKS_PER_SEC;

	cout << "The elapsed time during the program running: " << elapsed_secs<< endl;

	// clean the allocated memory
	if(weather){
		for(int i=0;i<height;i++){
			if(weather[i])
				delete [] weather[i];
		}
		delete [] weather;
	}

	return 0;
}
