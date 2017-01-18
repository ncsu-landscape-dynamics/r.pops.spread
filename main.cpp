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

#include <memory>
#include <stdexcept>

using namespace std;
using std::string;

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

Rtype radial_type_from_string(const string& text)
{
    if (text == "cauchy")
        return CAUCHY;
    else if (text == "cauchy_mix")
        return CAUCHY_MIX;
    else
        throw std::invalid_argument("Invalid value provided");
}

struct SodOptions
{
    struct Option *umca, *oaks, *lvtree, *ioaks;
    struct Option *nc_weather, *weather_value;
    struct Option *start_time, *end_time;
    struct Option *radial_type;
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
    opt.nc_weather->required = NO;

    opt.weather_value = G_define_option();
    opt.weather_value->type = TYPE_INTEGER;
    opt.weather_value->key = "weather_value";
    opt.weather_value->label = _("Value to be used as weather coeficient");
    opt.weather_value->description =
            _("Spatially and temporally constant weather coeficient"
              " (usually moisture times temperture in C)");
    opt.weather_value->required = NO;

    opt.start_time = G_define_option();
    opt.start_time->type = TYPE_INTEGER;
    opt.start_time->key = "start_time";
    opt.start_time->label = _("Start year for the simulation");
    opt.start_time->description = _("The first day of the year will be used");
    opt.start_time->required = YES;

    opt.end_time = G_define_option();
    opt.end_time->type = TYPE_INTEGER;
    opt.end_time->key = "end_time";
    opt.end_time->label = _("End year for the simulation");
    opt.end_time->description = _("The last day of the year will be used");
    opt.end_time->required = YES;

    opt.radial_type = G_define_option();
    opt.radial_type->type = TYPE_STRING;
    opt.radial_type->key = "radial_type";
    opt.radial_type->label = _("Radial distribution type");
    opt.radial_type->answer = "cauchy";
    opt.radial_type->options = "cauchy,cauchy_mix";

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

    // options for times are required ints
    int start_time = std::stoi(opt.start_time->answer);
    int end_time = std::stoi(opt.end_time->answer);
    if (start_time > end_time) {
        cerr << "Start date must precede the end date!!!" << endl;
        exit(EXIT_FAILURE);
    }


    std::shared_ptr<NcFile> weather_coeff = nullptr;
    double weather_value = 0;
    if (opt.nc_weather->answer)
        weather_coeff = std::make_shared<NcFile>(opt.nc_weather->answer, NcFile::ReadOnly);
    else if (opt.weather_value->answer)
        weather_value = std::stoi(opt.weather_value->answer);
    else
        weather_value = 1;  // no change (used in multiplication)

    if (weather_coeff && !weather_coeff->is_valid()) {
        cerr << "Can not open the weather coefficients file(.cn)!" << endl;
        exit(EXIT_FAILURE);
    }

    // TODO: do we need to free this? docs is 404
    NcVar *mcf_nc = nullptr;
    NcVar *ccf_nc = nullptr;

    if (weather_coeff && !(mcf_nc = weather_coeff->get_var("Mcoef"))) {
        cerr << "Can not read the moisture coefficients from the file!" <<
            endl;
        exit(EXIT_FAILURE);
    }

    if (weather_coeff && !(ccf_nc = weather_coeff->get_var("Ccoef"))) {
        cerr << "Can not read the temperature coefficients from the file!" <<
            endl;
        exit(EXIT_FAILURE);
    }

    double *mcf = nullptr;
    double *ccf = nullptr;
    double *weather = nullptr;
    if (weather_coeff) {
        mcf = new double[height * width];
        ccf = new double[height * width];
        weather = new double[height * width];
    }

    // Seasonality: Do you want the spread to be limited to certain months?
    bool ss = true;
    bool wind = true;

    // set the direction of the wind {N=0,NE=45,E=90,SE=135,S=180,SW=225,W=270,NW=315,NONE  // NONE means that there is no wind}
    Direction pwdir = NE;

    // set the spore rate
    double spore_rate = 4.4;
    Rtype rtype = radial_type_from_string(opt.radial_type->answer);
    double scale1 = 20.57;
    int kappa = 2;

    // initialize the start Date and end Date object
    Date dd_start(start_time, 01, 01);
    Date dd_end(end_time, 12, 31);

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

        if (weather_coeff) {
            // read the weather information
            if (!mcf_nc->set_cur(i,0,0)) {
                cerr << "Can not read the coefficients from the mcf_nc pointer to mcf array " << i << endl;
                exit(EXIT_FAILURE);
            }
            if (!ccf_nc->set_cur(i,0,0)) {
                cerr << "Can not read the coefficients from the ccf_nc pointer to ccf array "<< i << endl;
                exit(EXIT_FAILURE);
            }
            if (!mcf_nc->get(mcf, 1, height, width)) {
                cerr << "Can not get the record from mcf_nc " << i << endl;
                exit(EXIT_FAILURE);
            }
            if (!ccf_nc->get(ccf, 1, height, width)) {
                cerr << "Can not get the record from ccf_nc " << i << endl;
                exit(EXIT_FAILURE);
            }
            for (int j = 0; j < height; j++) {
                for (int k = 0; k < width; k++) {
                    weather[j * width + k] = mcf[j * width + k] * ccf[j * width + k];
                }
            }
        }

        // build the Sporulation object
        Sporulation sp1;
        sp1.SporeGen(I_umca_rast, weather, weather_value, spore_rate);

        sp1.SporeSpreadDisp(S_umca_rast, S_oaks_rast, I_umca_rast,
                            I_oaks_rast, lvtree_rast, rtype, weather,
                            weather_value, scale1, kappa, pwdir);

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
        delete[] mcf;
        delete[] ccf;
        delete[] weather;
    }

    return 0;
}
