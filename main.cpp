/*
 *
 *  SOD-model-cpp
 *
 *  Created on: Oct,2015
 *  Author: Zexi Chen(zchen22@ncsu.edu)
 *
 *
 */


#include <iostream>
#include "Img.h"
#include <netcdfcpp.h>

extern "C" {
#include <grass/gis.h>
#include <grass/glocale.h>
}

// to compile with netcdfcpp.h, add -lnetcdf_c++
#include "Spore.h"
//g++ -std=c++11 main.cpp Img.h Img.cpp Spore.h Spore.cpp -lgdal -lnetcdf_c++


using namespace std;
using std::string;

#define START_TIME 2000
#define END_TIME 2010
#define DIM 1

// Initialize infected trees for each species(!!Needed unless empirical info is available)
static Img initialize(Img& img1,Img& img2) {
    int re_width=0;
    int re_height=0;
    int **re_data=NULL;

    if (img1.getWidth() != img2.getWidth() ||
            img2.getHeight() != img2.getHeight()) {
        cerr << "The height or width of one image do not match with that of the other one!" << endl;
        return Img();
    }
    else {
        re_width = img1.getWidth();
        re_height = img1.getHeight();
        re_data = (int **)std::malloc(sizeof(int *) * re_height);
        int *stream = (int *)std::malloc(sizeof(int) * re_width * re_height);

        for (int i = 0; i < re_height; i++) {
            re_data[i] = &stream[i * re_width];
        }

        for (int i = 0; i < re_height; i++) {
            for (int j = 0; j < re_width; j++) {
                if (img2.data[i][j] > 0) {
                    if (img1.data[i][j] > img2.data[i][j])
                        re_data[i][j] =
                            img1.data[i][j] <
                            (img2.data[i][j] *
                             2) ? img1.data[i][j] : (img2.data[i][j] * 2);
                    else
                        re_data[i][j] = img1.data[i][j];
                }
                else {
                    re_data[i][j] = 0;
                }
            }
        }
        return Img(re_width, re_height, img1.getWEResolution(),
                   img1.getNSResolution(), re_data);
    }
}


// inputFname is used to retrieve GDAL information from the known input file
static void writeGeotiff(const char *inputFname, const char *outFname,
                         Img & img)
{
    // obtain information for output Geotiff images
    GDALDataset *inputDataset;
    GDALRasterBand *inputDataBand;

    GDALAllRegister();
    inputDataset = (GDALDataset *) GDALOpen(inputFname, GA_ReadOnly);
    double inputAdfGeoTransform[6];

    inputDataset->GetGeoTransform(inputAdfGeoTransform);

    // Setup driver
    const char *pszFormat = "GTiff";
    GDALDriver *gdalDriver;

    gdalDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    //char *pszSRS_WKT = NULL;
    int xSize = inputDataset->GetRasterXSize();
    int ySize = inputDataset->GetRasterYSize();
    GDALRasterBand *outBand;

    //Set output Dataset
    GDALDataset *outDataset;
    char **papszOptions = NULL;

    int *outstream = (int *)std::malloc(sizeof(int) * xSize * ySize);

    for (int i = 0; i < ySize; i++) {
        for (int j = 0; j < xSize; j++) {
            outstream[i * xSize + j] = img.data[i][j];
        }
    }

    // create output geotiff
    outDataset = gdalDriver->Create(outFname, xSize, ySize, 1, GDT_Byte, papszOptions);
    outDataset->SetGeoTransform(inputAdfGeoTransform);
    outDataset->SetProjection(inputDataset->GetProjectionRef());
    outBand = outDataset->GetRasterBand(1);
    outBand->RasterIO(GF_Write, 0, 0, xSize, ySize,outstream,
                      xSize, ySize, GDT_Int32, 0, 0);
    GDALClose((GDALDatasetH) outDataset);
    GDALClose((GDALDatasetH) inputDataset);
    if (outstream) {
        delete [] outstream;
    }
    CSLDestroy(papszOptions);
}

string generate_name(const string& basename, const Date& date)
{
    // counting on year being 4 digits
    auto year = G_double_to_basename_format(date.getYear(), 4, 0);
    auto month = G_double_to_basename_format(date.getMonth(), 2, 0);
    auto day = G_double_to_basename_format(date.getDay(), 2, 0);
    auto sep = G_get_basename_separator();
    string name = basename + sep + year + "_" + month + "_" + day;
    return name;
}

struct SodOptions
{
    struct Option *umca, *oaks, *lvtree, *ioaks;
    struct Option *nc_weather;
    struct Option *output, *output_series;
};

struct SodFlags
{
    //struct Flag *generate_seed;
};


int main(int argc, char *argv[])
{
    SodOptions opt;
    SodFlags flg;

    G_gisinit(argv[0]);

    struct GModule *module = G_define_module();

    G_add_keyword(_("raster"));
    G_add_keyword(_("spread"));
    G_add_keyword(_("model"));
    G_add_keyword(_("disease"));

    opt.umca = G_define_standard_option(G_OPT_R_INPUT);
    opt.umca->key = "umca";

    opt.oaks = G_define_standard_option(G_OPT_R_INPUT);
    opt.oaks->key = "oaks";

    opt.lvtree = G_define_standard_option(G_OPT_R_INPUT);
    opt.lvtree->key = "lvtree";

    opt.ioaks = G_define_standard_option(G_OPT_R_INPUT);
    opt.ioaks->key = "ioaks";

    opt.nc_weather = G_define_standard_option(G_OPT_F_INPUT);
    opt.nc_weather->key = "ncdf_weather";

    opt.output = G_define_standard_option(G_OPT_R_OUTPUT);

    opt.output_series = G_define_standard_option(G_OPT_R_OUTPUT);
    opt.output_series->key = "output_series";
    opt.output_series->required = NO;

    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

    // set the start point of the program
    clock_t begin = clock();

    // read the suspectible UMCA raster image
    Img umca_rast = Img::fromGrassRaster(opt.umca->answer);

    // read the SOD-affected oaks raster image
    Img oaks_rast = Img::fromGrassRaster(opt.oaks->answer);

    // read the living trees raster image
    Img lvtree_rast = Img::fromGrassRaster(opt.lvtree->answer);

    // read the initial infected oaks image
    Img I_oaks_rast = Img::fromGrassRaster(opt.ioaks->answer);

    // create the initial suspectible oaks image
    Img S_oaks_rast = oaks_rast - I_oaks_rast;

    // create the initial infected umca image
    Img I_umca_rast = initialize(umca_rast, I_oaks_rast);

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

    // retrieve the width and height of the images
    int width = umca_rast.getWidth();
    int height = umca_rast.getHeight();

    if (START_TIME > END_TIME) {
        cerr << "Start date must precede the end date!!!" << endl;
        exit(EXIT_FAILURE);
    }

    NcFile weatherCoeff(opt.nc_weather->answer, NcFile::ReadOnly);

    if (!weatherCoeff.is_valid()) {
        cerr << "Can not open the weather coefficients file(.cn)!" << endl;
        exit(EXIT_FAILURE);
    }

    NcVar *mcf_nc, *ccf_nc;

    if (!(mcf_nc = weatherCoeff.get_var("Mcoef"))) {
        cerr << "Can not read the moisture coefficients from the file!" <<
            endl;
        exit(EXIT_FAILURE);
    }

    if (!(ccf_nc = weatherCoeff.get_var("Ccoef"))) {
        cerr << "Can not read the temperature coefficients from the file!" <<
            endl;
        exit(EXIT_FAILURE);
    }

    double mcf[height][width];
    double ccf[height][width];
    double *weather = new double[height * width];

    // Seasonality: Do you want the spread to be limited to certain months?
    bool ss = true;
    bool wind = true;

    // set the direction of the wind {N=0,NE=45,E=90,SE=135,S=180,SW=225,W=270,NW=315,NONE  // NONE means that there is no wind}
    Direction pwdir = NE;

    // set the spore rate
    double spore_rate = 4.4;
    Rtype rtype = CAUCHY;
    double scale1 = 20.57;
    int kappa = 2;

    // initialize the start Date and end Date object
    Date dd_start(START_TIME, 01, 01);
    Date dd_end(END_TIME, 12, 31);

    // the variablbs created for the output to Geotiff file
    std::string s_year;
    std::string s_month;
    std::string s_day;
    // main simulation loop(weekly steps)
    for (int i = 0; dd_start.compareDate(dd_end);
         i++, dd_start.increasedByWeek()) {

        bool allInfected = true;

        for (int j = 0; j < height; j++) {
            for (int k = 0; k < width; k++) {
                if (S_oaks_rast.data[j][k] > 0)
                    allInfected = false;
            }
        }

        // if all the oaks are infected, then exit
        if (allInfected) {
            cerr << "In the " << dd_start.getYear() << "-" << dd_start.
                getMonth() << "-" << dd_start.
                getDay() << " all suspectible oaks are infected!" << endl;
            break;
        }

        // check whether the spore occurs in the month
        if (ss && dd_start.getMonth() > 9) {
            if (opt.output_series->answer) {
                // TODO: use end instead?
                string name = generate_name(opt.output_series->answer, dd_start);
                I_oaks_rast.toGrassRaster(name.c_str());
            }
            continue;
        }

        // read the weather information
        if (!mcf_nc->set_cur(i,0,0)) {
            cerr << "Can not read the coefficients from the mcf_nc pointer to mcf array " << i << endl;
            exit(EXIT_FAILURE);
        }

        if (!ccf_nc->set_cur(i,0,0)) {
            cerr << "Can not read the coefficients from the ccf_nc pointer to ccf array "<< i << endl;
            exit(EXIT_FAILURE);
        }

        if (!mcf_nc->get(&mcf[0][0], 1, height, width)) {
            cerr << "Can not get the record from mcf_nc " << i << endl;
            exit(EXIT_FAILURE);
        }

        if (!ccf_nc->get(&ccf[0][0], 1, height, width)) {
            cerr << "Can not get the record from ccf_nc " << i << endl;
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < height; j++) {
            for (int k = 0; k < width; k++) {
                weather[j * width + k] = mcf[j][k] * ccf[j][k];
            }
        }

        // build the Sporulation object
        Sporulation sp1;
        sp1.SporeGen(I_umca_rast,weather,spore_rate);

        sp1.SporeSpreadDisp(S_umca_rast, S_oaks_rast, I_umca_rast, I_oaks_rast, lvtree_rast,
                            rtype, weather, scale1, kappa, pwdir);

        s_year = std::to_string(dd_start.getYear());
        s_month = std::to_string(dd_start.getMonth());
        s_day = std::to_string(dd_start.getDay());

        if (opt.output_series->answer) {
            // TODO: use end instead?
            string name = generate_name(opt.output_series->answer, dd_start);
            I_oaks_rast.toGrassRaster(name.c_str());
        }
    }

    // write final result
    I_oaks_rast.toGrassRaster(opt.output->answer);

    // compute the time used when running the model
    clock_t end = clock();
    double elapsed_secs = double(end-begin) / CLOCKS_PER_SEC;

    cout << "The elapsed time during the program running: " << elapsed_secs << endl;

    // clean the allocated memory
    if (weather) {
        delete[] weather;
    }

    return 0;
}
