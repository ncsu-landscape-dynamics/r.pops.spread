/*
 *
 *  SOD-model-cpp
 *
 *  Created on: Oct,2015
 *  Author: Zexi Chen(zchen22@ncsu.edu)
 *
 *
 */


#include "Spore.h"
#include "Img.h"
#include "date.h"
#include "tcp_client.h"

extern "C" {
#include <grass/gis.h>
#include <grass/glocale.h>
}

#include <netcdfcpp.h>
#include <gdal/gdal.h>
#include <gdal/gdal_priv.h>

#include <atomic>
#include <arpa/inet.h> //inet_addr
#include <thread>

#include <map>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::atomic;
using std::thread;
using std::ref;


#define DIM 1

// Initialize infected trees for each species
// needed unless empirical info is available
static Img initialize(Img& img1,Img& img2) {
    if (img1.getWidth() != img2.getWidth() ||
            img2.getHeight() != img2.getHeight()) {
        cerr << "The height or width of one image do not match with that of the other one!" << endl;
        return Img();
    }
    else {
        auto re_width = img1.getWidth();
        auto re_height = img1.getHeight();
        auto out = Img(re_width, re_height, img1.getWEResolution(), img1.getNSResolution());

        for (int i = 0; i < re_height; i++) {
            for (int j = 0; j < re_width; j++) {
                if (img2(i, j) > 0) {
                    if (img1(i, j) > img2(i, j))
                        out(i, j) =
                            img1(i, j) <
                            (img2(i, j) *
                             2) ? img1(i, j) : (img2(i, j) * 2);
                    else
                        out(i, j) = img1(i, j);
                }
                else {
                    out(i, j) = 0;
                }
            }
        }
        return out;
    }
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

Direction direction_enum_from_string(const string& text)
{
    std::map<string, Direction> mapping{
        {"N", N}, {"NE", NE}, {"E", E}, {"SE", SE}, {"S", S},
        {"SW", SW}, {"W", W}, {"NW", NW}, {"NONE", NONE}
    };
    try {
        return mapping.at(text);
    }
    catch (const std::out_of_range&) {
        throw std::invalid_argument("direction_enum_from_string: Invalid"
                                    " value '" + text +"' provided");
    }
}

Rtype radial_type_from_string(const string& text)
{
    if (text == "cauchy")
        return CAUCHY;
    else if (text == "cauchy_mix")
        return CAUCHY_MIX;
    else
        throw std::invalid_argument("radial_type_from_string: Invalid"
                                    " value '" + text +"' provided");
}

bool seasonality_from_string(const string& text)
{
    if (text == "yes")
        return true;
    else if (text == "no")
        return false;
    else
        throw std::invalid_argument("seasonality_from_string: Invalid"
                                    " value '" + text +"' provided");
}

std::vector<double> weather_file_to_list(const string& filename)
{
    std::ifstream input(filename);
    std::vector<double> output;
    string line;
    while (std::getline(input, line))
    {
        double m, c;
        std::istringstream stream(line);
        stream >> m >> c;
        output.push_back(m * c);
    }
    return output;
}

bool all_infected(Img& S_oaks_rast)
{
    bool allInfected = true;
    for (int j = 0; j < S_oaks_rast.getHeight(); j++) {
        for (int k = 0; k < S_oaks_rast.getWidth(); k++) {
            if (S_oaks_rast(j, k) > 0)
                allInfected = false;
        }
    }
    return allInfected;
}

void get_spatial_weather(NcVar *mcf_nc, NcVar *ccf_nc, double* mcf, double* ccf, double* weather, int width, int height, int step)
{
    // read the weather information
    if (!mcf_nc->set_cur(step, 0, 0)) {
        cerr << "Can not read the coefficients from the mcf_nc pointer to mcf array " << step << endl;
        exit(EXIT_FAILURE);
    }
    if (!ccf_nc->set_cur(step, 0, 0)) {
        cerr << "Can not read the coefficients from the ccf_nc pointer to ccf array "<< step << endl;
        exit(EXIT_FAILURE);
    }
    if (!mcf_nc->get(mcf, 1, height, width)) {
        cerr << "Can not get the record from mcf_nc " << step << endl;
        exit(EXIT_FAILURE);
    }
    if (!ccf_nc->get(ccf, 1, height, width)) {
        cerr << "Can not get the record from ccf_nc " << step << endl;
        exit(EXIT_FAILURE);
    }
    for (int j = 0; j < height; j++) {
        for (int k = 0; k < width; k++) {
            weather[j * width + k] = mcf[j * width + k] * ccf[j * width + k];
        }
    }
}

void reload_UMCA_input(Img &umca, string map, Img &I_umca, Img &S_umca)
{
    umca = Img::fromGrassRaster(map.c_str());
    for (int i = 0; i < umca.getHeight(); i++) {
        for (int j = 0; j < umca.getWidth(); j++) {
            if (umca(i, j) < I_umca(i, j)) {
                I_umca(i, j) = umca(i, j);
            }
        }
    }
    S_umca = umca - I_umca;
}

void steering_client(tcp_client &c, string ip_address, int port, atomic<int> &instr_code, string &load_name, int &jump_date)
{
    int rec_error;
    string received;
    
    //connect to host
    c.conn(ip_address, port);

    while (true) {
        cout << "to receive" << endl;
        received = c.receive(200, rec_error);
        if (rec_error < 0){
            cerr << "receive failed\n";
            c.close_socket();
            instr_code.store(5);
            break;
        }
        else {
            if (received.substr(0, 3) == "cmd") {
                string cmd = received.substr(4, received.length() - 4);
                if (cmd == "play") {
                    instr_code.store(1);
                    cout << "play from " << endl;
                }
                else if (cmd == "pause") {
                    instr_code.store(2);
                    cout << "pause" << endl;
                }
                else if (cmd == "step") {
                    instr_code.store(3);
                    cout << "step" << endl;
                }
                else if (cmd == "stop") {
                    c.close_socket();
                    instr_code.store(5);
                    break;
                }
            } else if (received.substr(0, 4) == "load") {
                string name = received.substr(5, received.length() - 5);
                load_name = name;
                instr_code.store(4);
            } else if (received.substr(0, 4) == "move") {
                jump_date = std::stoi(received.substr(5, received.length() - 5));
            } else
                cout << "X" << received << "X" << rec_error << endl;
        }
    }
}

struct SodOptions
{
    struct Option *umca, *oaks, *lvtree, *ioaks;
    struct Option *nc_weather, *weather_value, *weather_file;
    struct Option *start_time, *end_time, *seasonality;
    struct Option *spore_rate, *wind;
    struct Option *radial_type, *scale_1, *scale_2, *kappa, *gamma;
    struct Option *seed;
    struct Option *output_series;
    struct Option *ip_address, *port;
};

struct SodFlags
{
    struct Flag *generate_seed;
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
    module->description = _("Stochastic landscape spread model of forest pathogen - Sudden Oak Death (SOD)");

    opt.umca = G_define_standard_option(G_OPT_R_INPUT);
    opt.umca->key = "umca";
    opt.umca->description = _("Input bay laurel (UMCA) raster map");
    opt.umca->guisection = _("Input");

    opt.oaks = G_define_standard_option(G_OPT_R_INPUT);
    opt.oaks->key = "oaks";
    opt.oaks->description = _("Input SOD-oaks raster map");
    opt.oaks->guisection = _("Input");

    opt.lvtree = G_define_standard_option(G_OPT_R_INPUT);
    opt.lvtree->key = "lvtree";
    opt.lvtree->description = _("Input live tree (all) raster map");
    opt.lvtree->guisection = _("Input");

    // TODO: is this oaks?
    opt.ioaks = G_define_standard_option(G_OPT_R_INPUT);
    opt.ioaks->key = "ioaks";
    opt.ioaks->description = _("Initial sources of infection raster map");
    opt.ioaks->guisection = _("Input");

    opt.output_series = G_define_standard_option(G_OPT_R_BASENAME_OUTPUT);
    opt.output_series->key = "output_series";
    opt.output_series->description = _("Basename for output series");
    opt.output_series->required = NO;
    opt.output_series->guisection = _("Output");

    opt.wind = G_define_option();
    opt.wind->type = TYPE_STRING;
    opt.wind->key = "wind";
    opt.wind->label = _("Prevailing wind direction");
    opt.wind->description = _("NONE means that there is no wind");
    opt.wind->options = "N,NE,E,SE,S,SW,W,NW,NONE";
    opt.wind->required = YES;
    opt.wind->guisection = _("Weather");

    opt.nc_weather = G_define_standard_option(G_OPT_F_BIN_INPUT);
    opt.nc_weather->key = "ncdf_weather";
    opt.nc_weather->description = _("Weather data");
    opt.nc_weather->required = NO;
    opt.nc_weather->guisection = _("Weather");

    opt.weather_file = G_define_standard_option(G_OPT_F_INPUT);
    opt.weather_file->key = "weather_file";
    opt.weather_file->label = _("Text file with weather");
    opt.weather_file->description =
            _("Moisture and temperature");
    opt.weather_file->required = NO;
    opt.weather_file->guisection = _("Weather");

    opt.weather_value = G_define_option();
    opt.weather_value->type = TYPE_INTEGER;
    opt.weather_value->key = "weather_value";
    opt.weather_value->label = _("Value to be used as weather coeficient");
    opt.weather_value->description =
            _("Spatially and temporally constant weather coeficient"
              " (usually moisture times temperture)");
    opt.weather_value->required = NO;
    opt.weather_value->guisection = _("Weather");

    opt.start_time = G_define_option();
    opt.start_time->type = TYPE_INTEGER;
    opt.start_time->key = "start_time";
    opt.start_time->label = _("Start year for the simulation");
    opt.start_time->description = _("The first day of the year will be used");
    opt.start_time->required = YES;
    opt.start_time->guisection = _("Time");

    opt.end_time = G_define_option();
    opt.end_time->type = TYPE_INTEGER;
    opt.end_time->key = "end_time";
    opt.end_time->label = _("End year for the simulation");
    opt.end_time->description = _("The last day of the year will be used");
    opt.end_time->required = YES;
    opt.end_time->guisection = _("Time");

    opt.seasonality = G_define_option();
    opt.seasonality->type = TYPE_STRING;
    opt.seasonality->key = "seasonality";
    opt.seasonality->label = _("Seasonal spread");
    opt.seasonality->description = _("Spread limited to certain months (season)");
    opt.seasonality->options = "yes,no";
    opt.seasonality->answer = "yes";
    opt.seasonality->guisection = _("Time");

    opt.spore_rate = G_define_option();
    opt.spore_rate->type = TYPE_DOUBLE;
    opt.spore_rate->key = "spore_rate";
    opt.spore_rate->label = _("Spore production rate per week for each infected tree");
    opt.spore_rate->answer = "4.4";
    opt.spore_rate->guisection = _("Spores");

    opt.radial_type = G_define_option();
    opt.radial_type->type = TYPE_STRING;
    opt.radial_type->key = "radial_type";
    opt.radial_type->label = _("Radial distribution type");
    opt.radial_type->answer = "cauchy";
    opt.radial_type->options = "cauchy,cauchy_mix";
    opt.radial_type->guisection = _("Spores");

    opt.scale_1 = G_define_option();
    opt.scale_1->type = TYPE_DOUBLE;
    opt.scale_1->key = "scale_1";
    opt.scale_1->label = _("Scale parameter for the first Cauchy distribution");
    opt.scale_1->answer = "20.57";
    opt.scale_1->guisection = _("Spores");

    opt.scale_2 = G_define_option();
    opt.scale_2->type = TYPE_DOUBLE;
    opt.scale_2->key = "scale_2";
    opt.scale_2->label = _("Scale parameter for the second Cauchy distribution");
    opt.scale_2->guisection = _("Spores");

    opt.kappa = G_define_option();
    opt.kappa->type = TYPE_DOUBLE;
    opt.kappa->key = "kappa";
    opt.kappa->label = _("Concentration parameter for the von Mises distribution");
    opt.kappa->answer = "2";
    opt.kappa->guisection = _("Spores");

    opt.gamma = G_define_option();
    opt.gamma->type = TYPE_DOUBLE;
    opt.gamma->key = "gamma";
    opt.gamma->label = _("Gamma parameter for Bernoulli distribution");
    opt.gamma->description = _("Probability of using the first Cauchy distribution");
    opt.gamma->options = "0-1";
    opt.gamma->guisection = _("Spores");

    opt.seed = G_define_option();
    opt.seed->key = "random_seed";
    opt.seed->type = TYPE_INTEGER;
    opt.seed->required = NO;
    opt.seed->label = _("Seed for random number generator");
    opt.seed->description =
        _("The same seed can be used to obtain same results"
          " or random seed can be generated by other means.");
    opt.seed->guisection = _("Randomness");

    flg.generate_seed = G_define_flag();
    flg.generate_seed->key = 's';
    flg.generate_seed->label =
        _("Generate random seed (result is non-deterministic)");
    flg.generate_seed->description =
        _("Automatically generates random seed for random number"
          " generator (use when you don't want to provide the seed option)");
    flg.generate_seed->guisection = _("Randomness");
    
    opt.ip_address = G_define_option();
    opt.ip_address->key = "ip_address";
    opt.ip_address->type = TYPE_STRING;
    opt.ip_address->required = YES;
    opt.ip_address->description = _("IP address of steering server");
    opt.ip_address->guisection = _("Steering");

    opt.port = G_define_option();
    opt.port->key = "port";
    opt.port->type = TYPE_INTEGER;
    opt.port->required = YES;
    opt.port->description = _("Port of steering server");
    opt.port->guisection = _("Steering");

    G_option_exclusive(opt.seed, flg.generate_seed, NULL);
    G_option_required(opt.seed, flg.generate_seed, NULL);

    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

    // set the start point of the program
    clock_t begin = clock();

    // Seasonality: Do you want the spread to be limited to certain months?
    bool ss = seasonality_from_string(opt.seasonality->answer);

    Direction pwdir = direction_enum_from_string(opt.wind->answer);

    // set the spore rate
    double spore_rate = std::stod(opt.spore_rate->answer);
    Rtype rtype = radial_type_from_string(opt.radial_type->answer);
    double scale1 = std::stod(opt.scale_1->answer);
    double scale2 = 0;
    if (rtype == CAUCHY_MIX && !opt.scale_2->answer)
        G_fatal_error(_("The option %s is required for %s=%s"),
                      opt.scale_2->key, opt.radial_type->key,
                      opt.radial_type->answer);
    else if (opt.scale_2->answer)
        scale2 = std::stod(opt.scale_2->answer);
    double kappa = std::stod(opt.kappa->answer);
    double gamma = 0.0;
    if (rtype == CAUCHY_MIX && !opt.gamma->answer)
        G_fatal_error(_("The option %s is required for %s=%s"),
                      opt.gamma->key, opt.radial_type->key,
                      opt.radial_type->answer);
    else if (opt.gamma->answer)
        gamma = std::stod(opt.gamma->answer);

    // initialize the start Date and end Date object
    // options for times are required ints
    int start_time = std::stoi(opt.start_time->answer);
    int end_time = std::stoi(opt.end_time->answer);
    if (start_time > end_time) {
        cerr << "Start date must precede the end date!!!" << endl;
        exit(EXIT_FAILURE);
    }
    Date dd_start(start_time, 01, 01);
    Date dd_end(end_time, 12, 31);

    unsigned seed_value;
    if (opt.seed->answer) {
        seed_value = std::stoul(opt.seed->answer);
        G_verbose_message(_("Read random seed from %s option: %ud"),
                          opt.seed->key, seed_value);
    } else {
        // flag of option is required
        std::random_device rd;
        seed_value = rd();
        G_verbose_message(_("Generated random seed (-%c): %ud"),
                          flg.generate_seed->key, seed_value);
    }

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

    // retrieve the width and height of the images
    int width = umca_rast.getWidth();
    int height = umca_rast.getHeight();

    std::shared_ptr<NcFile> weather_coeff = nullptr;
    std::vector<double> weather_values;
    double weather_value = 0;
    if (opt.nc_weather->answer)
        weather_coeff = std::make_shared<NcFile>(opt.nc_weather->answer, NcFile::ReadOnly);
    else if (opt.weather_file->answer)
        weather_values = weather_file_to_list(opt.weather_file->answer);
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

    // setup steering variables
    Date dd_current(dd_start);
    Date dd_current_end(dd_start);

    // setup client
    tcp_client c;
    atomic<int> instr_code(0);
    thread client_thread;
    string ip = string(opt.ip_address->answer);
    string load_name = "";
    int jump_date = 0;
    int port = atoi(opt.port->answer);
    client_thread = thread(steering_client, ref(c), ip, port, ref(instr_code), ref(load_name), ref(jump_date));

    // build the Sporulation object
    Sporulation sp1(seed_value, I_umca_rast);

    // main simulation loop (weekly steps)
    int i = 0;
    while (true) {
        int code = instr_code.load();
        instr_code.store(0);
        if (code == 1) { // play from
            dd_current_end = dd_end;
        }
        else if (code == 2) { // pause
            dd_current_end = dd_current.getYearEnd();
        }
        else if (code == 3) { // 1 step forward
            dd_current_end = dd_current.getNextYearEnd();
            if (dd_current_end > dd_end)
                dd_current_end = dd_end;
        }
        else if (code == 4) { // load data
            cout << "loading data: " << load_name << endl;
            reload_UMCA_input(umca_rast, load_name, I_umca_rast, S_umca_rast);
        }
        else if (code == 5) { // complete stop
            break;
        }
        
        if (dd_current_end > dd_start && dd_current <= dd_current_end) {
            cout << i << endl;
            // here do spread
            // if all the oaks are infected, then exit
            if (all_infected(S_oaks_rast)) {
                cerr << "In the " << dd_start.getYear() << "-" << dd_start.
                        getMonth() << "-" << dd_start.
                        getDay() << " all suspectible oaks are infected!" << endl;
                continue;
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
                get_spatial_weather(mcf_nc, ccf_nc, mcf, ccf, weather, width, height, i);
            } else if (!weather_values.empty()) {
                weather_value = weather_values[i];
            }

            sp1.SporeGen(I_umca_rast, weather, weather_value, spore_rate);

            sp1.SporeSpreadDisp(S_umca_rast, S_oaks_rast, I_umca_rast,
                                I_oaks_rast, lvtree_rast, rtype, weather,
                                weather_value, scale1, kappa, pwdir, scale2,
                                gamma);

            if (dd_current.isYearEnd()) {
                if (opt.output_series->answer) {
                    string name = generate_name(opt.output_series->answer, dd_current);
                    I_oaks_rast.toGrassRaster(name.c_str());
                }
                c.send_data("output");
            }

            dd_current.increasedByWeek();
            i += 1;
        }
        else {
            // paused
        }
    }

    client_thread.join();
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
