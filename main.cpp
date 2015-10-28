#include <iostream>
#include "Img.h"
#include <netcdfcpp.h>
#include "Spore.h"
//#include <opencv2/highgui/highgui.hpp>
// to compile, add g++ main.cpp -ldgal
//#include <gdal/gdal_priv.h>
//http://www.unidata.ucar.edu/software/netcdf/examples/programs/
//http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial/
//<netcdfcpp.h>
// to compile with netcdfcpp.h, add -lnetcdf_c++
// http://docs.opencv.org/modules/core/doc/intro.html
//g++ -std=c++11 main.cpp Img.h Img.cpp Spore.h Spore.cpp -lgdal -lnetcdf_c++


using namespace std;

#define UMCA_RAST_PATH "./layers/UMCA_den_100m.img"
#define OAKS_RAST_PATH "./layers/OAKS_den_100m.img"
#define LVTREE_RAST_PATH "./layers/TPH_den_100m.img"
#define I_OAKS_RAST_PATH "./layers/init_2000_cnt.img"
#define BACKGROUND_IMAGE_PATH "./layers/ortho_5m_color.tif"
#define START_TIME 2000
#define END_TIME 2007
#define WEATHER_COEFF_PATH "./layers/weather/weatherCoeff_2000_2007.nc"
#define DIM 1

int main()
{
	clock_t begin = clock();
	// read the UMCA raster image
	Img umca_rast(UMCA_RAST_PATH);
	//umca_rast.readImage(UMCA_RAST_PATH);
	// read the oaks raster image

	Img oaks_rast(OAKS_RAST_PATH);
	//oaks_rast.readImage(OAKS_RAST_PATH);
	// read the living trees raster image
	Img lvtree_rast(LVTREE_RAST_PATH);
	//lvtree_rast.readImage(LVTREE_RAST_PATH);

	// SOD-immune trees image
	Img SOD_rast = umca_rast + oaks_rast;
	Img IMM_rast = lvtree_rast - SOD_rast;

	Img I_oaks_rast(I_OAKS_RAST_PATH);

	Img I_umca_rast = I_oaks_rast * 2;

	Img S_oaks_rast = oaks_rast - I_oaks_rast;
	Img S_umca_rast = umca_rast - I_umca_rast;
/*
	for(int i=0;i<S_oaks_rast.getHeight();i++){
		for(int j=0;j<S_oaks_rast.getWidth();j++){
			cout << S_oaks_rast.data[i][j] << " ";
		}
		cout << endl;
	}
*/
	Img bkr_img(BACKGROUND_IMAGE_PATH);

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

	bool wind = true;
	Direction pwdir=NE;
	double spore_rate = 4.4;
	int nth_output = 4;
	Rtype rtype = CAUCHY;
	double scale1 = 20.57;
	int kappa = 2;
	for(int i=0;i<=415;i++){
		bool allInfected = true;
		for(int j=0;j<height;j++){
			for(int k=0;k<width;k++){
				if(S_oaks_rast.data[j][k]>0)
					allInfected = false;
			}
		}

		if(allInfected){
			cerr << "In the " << i << "th iteration, all suspectible oaks are infected!" << endl;
			break;
		}
	// we can set the first dimension of the set_cur(first,second,third), which is the timestamp
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

		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				weather[i][j] = mcf[i][j] * ccf[i][j];
			}
		}

		Sporulation sp1;
		sp1.SporeGen(I_umca_rast,weather,spore_rate);
		sp1.SporeSpreadDisp(S_umca_rast, S_oaks_rast, I_umca_rast, I_oaks_rast, IMM_rast, 
				rtype, weather, scale1, kappa, pwdir);

		cout << "----------"  << i << "-------------" << endl;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				cout << I_oaks_rast.data[i][j] << " ";
			}
			cout << endl;
		}
	}

	clock_t end = clock();
	double elapsed_secs = double(end-begin) / CLOCKS_PER_SEC;

	cout << "The elapsed time during the program running: " << elapsed_secs<< endl;

	/*
	for(int i=0;i<bkr_img.getHeight();i++){
		for(int j=0;j<bkr_img.getWidth();j++){
			cout << bkr_img.data[i*bkr_img.getWidth() + j]<< " ";
		}
		cout << endl;
	}
	*/

	if(weather){
		for(int i=0;i<height;i++){
			if(weather[i])
				delete [] weather[i];
		}
		delete [] weather;
	}

	return 0;
}

	/*
	GDALDataset *dataset;
	GDALRasterBand *dataBand;
	int width;
	int height;
	int *data;
	GDALAllRegister();
	dataset = (GDALDataset *)GDALOpen(PATH, GA_ReadOnly);
	if(!dataset ){
		cerr << "Can not read the image!" << endl;
	}

	width = dataset->GetRasterXSize();
	height = dataset->GetRasterYSize();

	dataBand = dataset->GetRasterBand(1);
	data = (int *) CPLMalloc(sizeof(int)*width*height);
	dataBand->RasterIO(GF_Read,0,0,width,height,data,
		width,height,GDT_Int32,0,0);
	*/

/*
	for(int i=0;i<height;i++){
		for(int j=0;j<width;j++){
				cout << mcf[i][j] << " ";
		}
		cout << endl;
	}
*/