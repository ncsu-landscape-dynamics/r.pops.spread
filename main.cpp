/*
 * SOD model
 *
 * Copyright (C) 2015-2017 by the authors.
 *
 * Authors: Zexi Chen (zchen22 ncsu edu)
 *          Vaclav Petras (wenzeslaus gmail com)
 *          Anna Petrasova (kratochanna gmail com)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#include "graster.hpp"

#include "pops/model.hpp"
#include "pops/model_type.hpp"
#include "pops/date.hpp"
#include "pops/raster.hpp"
#include "pops/kernel.hpp"
#include "pops/treatments.hpp"
#include "pops/spread_rate.hpp"
#include "pops/statistics.hpp"
#include "pops/scheduling.hpp"
#include "pops/quarantine.hpp"
#include "pops/host_pool.hpp"
#include "pops/pest_pool.hpp"
#include "pops/multi_host_pool.hpp"

extern "C" {
#include <grass/gis.h>
#include <grass/glocale.h>
#include <grass/vector.h>
#include <grass/raster.h>
}

#include <map>
#include <tuple>
#include <vector>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include <sys/stat.h>

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::round;
using std::isnan;

using namespace pops;

#define DIM 1

void fatal_option_required_for_value(struct Option* required, struct Option* given)
{
    G_fatal_error(
        _("The option %s is required for %s=%s"),
        required->key,
        given->key,
        given->answer);
}

// check if a file exists
inline bool file_exists(const char* name)
{
    struct stat buffer;
    return (stat(name, &buffer) == 0);
}

inline void file_exists_or_fatal_error(struct Option* option)
{
    if (option->answer && !file_exists(option->answer))
        G_fatal_error(
            _("Option %s: File %s does not exist"), option->key, option->answer);
}

string generate_name(const string& basename, const Date& date)
{
    // counting on year being 4 digits
    auto year = G_double_to_basename_format(date.year(), 4, 0);
    auto month = G_double_to_basename_format(date.month(), 2, 0);
    auto day = G_double_to_basename_format(date.day(), 2, 0);
    auto sep = G_get_basename_separator();
    string name = basename + sep + year + "_" + month + "_" + day;
    return name;
}

void write_average_area(
    const std::vector<Img>& infected,
    const char* raster_name,
    double ew_res,
    double ns_res,
    const std::vector<std::vector<int>>& suitable_cells)
{
    struct History hist;
    double avg = 0;
    for (unsigned i = 0; i < infected.size(); i++) {
        avg += area_of_infected(infected[i], ew_res, ns_res, suitable_cells);
    }
    avg /= infected.size();
    string avg_string = "Average infected area: " + std::to_string(avg);
    Rast_read_history(raster_name, "", &hist);
    Rast_set_history(&hist, HIST_KEYWRD, avg_string.c_str());
    Rast_write_history(raster_name, &hist);
}

inline Date treatment_date_from_string(const string& text)
{
    try {
        return Date(text);
    }
    catch (std::invalid_argument&) {
        G_fatal_error(_("Date <%s> is invalid"), text.c_str());
    }
}

inline void seasonality_from_option(Config& config, const Option* opt)
{
    config.set_season_start_end_month(opt->answers[0], opt->answers[1]);
}

unsigned int get_num_answers(struct Option* opt)
{
    unsigned int i = 0;
    if (opt->answers)
        for (i = 0; opt->answers[i]; i++)
            ;
    return i;
}

/** Read list of names from a file
 *
 * Assumes one name per line. The current implementation basically splits the file by
 * lines and returns the resulting list of strings.
 */
std::vector<string> read_names(const char* filename)
{
    std::vector<string> names;
    std::ifstream file(filename);
    string line;
    while (std::getline(file, line)) {
        names.push_back(line);
    }
    return names;
}

/** From a file option representing a raster series, create series of rasters.
 *
 * Calls fatal error if the input is not correct. It will read all rasters if there is
 * more than needed.
 */
std::vector<DImg> file_to_series(struct Option* option, int expected)
{
    file_exists_or_fatal_error(option);
    auto names{read_names(option->answer)};
    if (names.size() < size_t(expected))
        G_fatal_error(
            _("Not enough names for %s in '%s' (%d expected)"),
            option->key,
            option->answer,
            expected);
    std::vector<DImg> rasters;
    for (const auto& name : names) {
        rasters.push_back(raster_from_grass_float(name));
    }
    return rasters;
}

/*!
 * Warns about depreciated option value
 *
 * It uses the answer member. If the answer is not set,
 * nothing is tested.
 *
 * \param opt Pointer to a valid option structure
 * \param depreciated Value which is depreciated
 * \param current Value which should be used instead
 */
void warn_about_depreciated_option_value(
    const Option* opt, const string& depreciated, const string& current)
{
    if (opt->answer && opt->answer == depreciated) {
        G_warning(
            _("The value <%s> for option %s is depreciated."
              " Use value <%s> instead."),
            opt->answer,
            opt->key,
            current.c_str());
    }
}

std::vector<double> weather_file_to_list(const string& filename)
{
    std::ifstream input(filename);
    std::vector<double> output;
    string line;
    while (std::getline(input, line)) {
        double m, c;
        std::istringstream stream(line);
        stream >> m >> c;
        output.push_back(m * c);
    }
    return output;
}

/** Checks if there are any susceptible hosts left */
bool all_infected(Img& susceptible)
{
    for (Img::IndexType j = 0; j < susceptible.rows(); j++)
        for (Img::IndexType k = 0; k < susceptible.cols(); k++)
            if (susceptible(j, k) > 0)
                return false;
    return true;
}

struct PoPSOptions
{
    struct Option *host, *total_plants, *infected, *outside_spores;
    struct Option* model_type;
    struct Option* latency_period;
    struct Option* dispersers_to_soils;
    struct Option* soil_survival_steps;
    struct Option *moisture_coefficient_file, *temperature_coefficient_file;
    struct Option* weather_coefficient_file;
    struct Option* weather_coefficient_stddev_file;
    struct Option* survival_rate_month;
    struct Option* survival_rate_day;
    struct Option* survival_rate_file;
    struct Option *lethal_temperature, *lethal_temperature_months;
    struct Option* temperature_file;
    struct Option *start_date, *end_date, *seasonality;
    struct Option *step_unit, *step_num_units;
    struct Option* treatments;
    struct Option *treatment_date, *treatment_length;
    struct Option* treatment_app;
    struct Option* reproductive_rate;
    struct Option *natural_kernel, *natural_scale;
    struct Option *natural_direction, *natural_kappa;
    struct Option *anthro_kernel, *anthro_scale;
    struct Option *anthro_direction, *anthro_kappa;
    struct Option* percent_natural_dispersal;
    struct Option *infected_to_dead_rate, *first_year_to_die;
    struct Option *mortality_frequency, *mortality_frequency_n;
    struct Option* dead_series;
    struct Option* seed;
    struct Option* seeds;
    struct Option* runs;
    struct Option* threads;
    struct Option* single_series;
    struct Option *average, *average_series;
    struct Option *stddev, *stddev_series;
    struct Option *probability, *probability_series;
    struct Option* spread_rate_output;
    struct Option *quarantine, *quarantine_output, *quarantine_directions;
    struct Option *dispersers_output, *established_dispersers_output;
    struct Option *output_frequency, *output_frequency_n;
};

struct PoPSFlags
{
    struct Flag* mortality;
    struct Flag* generate_seed;
};

int main(int argc, char* argv[])
{
    PoPSOptions opt;
    PoPSFlags flg;

    G_gisinit(argv[0]);

    struct GModule* module = G_define_module();

    G_add_keyword(_("raster"));
    G_add_keyword(_("spread"));
    G_add_keyword(_("model"));
    G_add_keyword(_("simulation"));
    G_add_keyword(_("disease"));
    G_add_keyword(_("pest"));
    module->description =
        _("A dynamic species distribution model for pest or "
          "pathogen spread in forest or agricultural ecosystems (PoPS)");

    opt.host = G_define_standard_option(G_OPT_R_INPUT);
    opt.host->key = "host";
    opt.host->label = _("Input host raster map");
    opt.host->description = _("Number of hosts per cell.");
    opt.host->guisection = _("Input");

    opt.total_plants = G_define_standard_option(G_OPT_R_INPUT);
    opt.total_plants->key = "total_plants";
    opt.total_plants->label = _("Input raster map of total plants");
    opt.total_plants->description = _("Number of all plants per cell");
    opt.total_plants->guisection = _("Input");

    opt.infected = G_define_standard_option(G_OPT_R_INPUT);
    opt.infected->key = "infected";
    opt.infected->label = _("Input raster map of initial infection");
    opt.infected->description = _("Number of infected hosts per cell");
    opt.infected->guisection = _("Input");

    opt.average = G_define_standard_option(G_OPT_R_OUTPUT);
    opt.average->key = "average";
    opt.average->description = _("Average infected across multiple runs");
    opt.average->guisection = _("Output");
    opt.average->required = NO;

    opt.average_series = G_define_standard_option(G_OPT_R_BASENAME_OUTPUT);
    opt.average_series->key = "average_series";
    opt.average_series->description =
        _("Basename for output series of average infected across multiple runs");
    opt.average_series->required = NO;
    opt.average_series->guisection = _("Output");

    opt.single_series = G_define_standard_option(G_OPT_R_BASENAME_OUTPUT);
    opt.single_series->key = "single_series";
    opt.single_series->description =
        _("Basename for output series of infected as single stochastic run");
    opt.single_series->required = NO;
    opt.single_series->guisection = _("Output");

    opt.stddev = G_define_standard_option(G_OPT_R_OUTPUT);
    opt.stddev->key = "stddev";
    opt.stddev->description = _("Standard deviations");
    opt.stddev->required = NO;
    opt.stddev->guisection = _("Output");

    opt.stddev_series = G_define_standard_option(G_OPT_R_BASENAME_OUTPUT);
    opt.stddev_series->key = "stddev_series";
    opt.stddev_series->description =
        _("Basename for output series of standard deviations");
    opt.stddev_series->required = NO;
    opt.stddev_series->guisection = _("Output");

    opt.probability = G_define_standard_option(G_OPT_R_OUTPUT);
    opt.probability->key = "probability";
    opt.probability->description = _("Infection probability (in percent)");
    opt.probability->required = NO;
    opt.probability->guisection = _("Output");

    opt.probability_series = G_define_standard_option(G_OPT_R_BASENAME_OUTPUT);
    opt.probability_series->key = "probability_series";
    opt.probability_series->description =
        _("Basename for output series of probabilities");
    opt.probability_series->required = NO;
    opt.probability_series->guisection = _("Output");

    opt.outside_spores = G_define_standard_option(G_OPT_V_OUTPUT);
    opt.outside_spores->key = "outside_spores";
    opt.outside_spores->description =
        _("Output vector map of spores or pest units outside of modeled area");
    opt.outside_spores->required = NO;
    opt.outside_spores->guisection = _("Output");

    opt.spread_rate_output = G_define_standard_option(G_OPT_F_OUTPUT);
    opt.spread_rate_output->key = "spread_rate_output";
    opt.spread_rate_output->description =
        _("Output CSV file containg yearly spread rate in N, S, E, W directions");
    opt.spread_rate_output->required = NO;
    opt.spread_rate_output->guisection = _("Output");

    opt.quarantine = G_define_standard_option(G_OPT_R_INPUT);
    opt.quarantine->key = "quarantine";
    opt.quarantine->required = NO;
    opt.quarantine->label = _("Input raster map of quarantine areas");
    opt.quarantine->description = _("Values > 0 represent quarantine areas");
    opt.quarantine->guisection = _("Input");

    opt.quarantine_output = G_define_standard_option(G_OPT_F_OUTPUT);
    opt.quarantine_output->key = "quarantine_output";
    opt.quarantine_output->description =
        _("Output CSV file containg yearly quarantine information");
    opt.quarantine_output->required = NO;
    opt.quarantine_output->guisection = _("Output");

    opt.quarantine_directions = G_define_option();
    opt.quarantine_directions->type = TYPE_STRING;
    opt.quarantine_directions->key = "quarantine_directions";
    opt.quarantine_directions->label = _("Quarantine directions to consider");
    opt.quarantine_directions->description =
        _("Comma separated directions to include"
          "in the quarantine direction analysis, e.g., 'N,E' "
          "(by default all directions (N, S, E, W) are considered)");
    opt.quarantine_directions->required = NO;
    opt.quarantine_directions->guisection = _("Output");

    opt.model_type = G_define_option();
    opt.model_type->type = TYPE_STRING;
    opt.model_type->key = "model_type";
    opt.model_type->label = _("Epidemiological model type");
    opt.model_type->answer = const_cast<char*>("SI");
    opt.model_type->options = "SI,SEI";
    opt.model_type->descriptions =
        _("SI;Susceptible-infected epidemiological model;"
          "SEI;Susceptible-exposed-infected epidemiological model"
          " (uses latency_period)");
    opt.model_type->required = YES;
    opt.model_type->guisection = _("Model");

    opt.latency_period = G_define_option();
    opt.latency_period->type = TYPE_INTEGER;
    opt.latency_period->key = "latency_period";
    opt.latency_period->label = _("Latency period in simulation steps");
    opt.latency_period->description =
        _("How long it takes for a hosts to become infected after being exposed"
          " (unit is a simulation step)");
    opt.latency_period->required = NO;
    opt.latency_period->guisection = _("Model");

    opt.dispersers_to_soils = G_define_option();
    opt.dispersers_to_soils->type = TYPE_DOUBLE;
    opt.dispersers_to_soils->key = "dispersers_to_soils";
    opt.dispersers_to_soils->label = _("Ratio of dispersers going into soil");
    opt.dispersers_to_soils->description =
        _("Ratio (percentage) of generated dispersers going into soil instead "
          "of being dispersed by the kernel");
    opt.dispersers_to_soils->options = "0-1";
    opt.dispersers_to_soils->required = NO;
    opt.dispersers_to_soils->guisection = _("Model");

    opt.soil_survival_steps = G_define_option();
    opt.soil_survival_steps->type = TYPE_INTEGER;
    opt.soil_survival_steps->key = "soil_survival_steps";
    opt.soil_survival_steps->label = _("Steps dispersers stay in the soil");
    opt.soil_survival_steps->description =
        _("Number of simulation steps dispersers survive in the soil");
    opt.soil_survival_steps->options = "1-";
    opt.soil_survival_steps->required = NO;
    opt.soil_survival_steps->guisection = _("Model");

    opt.treatments = G_define_standard_option(G_OPT_R_INPUT);
    opt.treatments->key = "treatments";
    opt.treatments->multiple = YES;
    opt.treatments->description =
        _("Raster map(s) of treatments (treated 1, otherwise 0)");
    opt.treatments->required = NO;
    opt.treatments->guisection = _("Treatments");

    opt.treatment_date = G_define_option();
    opt.treatment_date->key = "treatment_date";
    opt.treatment_date->type = TYPE_STRING;
    opt.treatment_date->multiple = YES;
    opt.treatment_date->description =
        _("Dates when treatments are applied (e.g. 2020-01-15)");
    opt.treatment_date->required = NO;
    opt.treatment_date->guisection = _("Treatments");

    opt.treatment_length = G_define_option();
    opt.treatment_length->type = TYPE_INTEGER;
    opt.treatment_length->key = "treatment_length";
    opt.treatment_length->multiple = YES;
    opt.treatment_length->label = _("Treatment length in days");
    opt.treatment_length->description =
        _("Treatment length 0 results in simple removal of host, length > 0 makes"
          " host resistant for certain number of days");
    opt.treatment_length->required = NO;
    opt.treatment_length->guisection = _("Treatments");

    opt.treatment_app = G_define_option();
    opt.treatment_app->key = "treatment_application";
    opt.treatment_app->type = TYPE_STRING;
    opt.treatment_app->multiple = NO;
    opt.treatment_app->description = _("Type of treatmet application");
    opt.treatment_app->options = "ratio_to_all,all_infected_in_cell";
    opt.treatment_app->required = NO;
    opt.treatment_app->answer = const_cast<char*>("ratio_to_all");
    opt.treatment_app->guisection = _("Treatments");

    opt.moisture_coefficient_file = G_define_standard_option(G_OPT_F_INPUT);
    opt.moisture_coefficient_file->key = "moisture_coefficient_file";
    opt.moisture_coefficient_file->label =
        _("Input file with one moisture coefficient map name per line");
    opt.moisture_coefficient_file->description = _("Moisture coefficient");
    opt.moisture_coefficient_file->required = NO;
    opt.moisture_coefficient_file->guisection = _("Weather");

    opt.temperature_coefficient_file = G_define_standard_option(G_OPT_F_INPUT);
    opt.temperature_coefficient_file->key = "temperature_coefficient_file";
    opt.temperature_coefficient_file->label =
        _("Input file with one temperature coefficient map name per line");
    opt.temperature_coefficient_file->description = _("Temperature coefficient");
    opt.temperature_coefficient_file->required = NO;
    opt.temperature_coefficient_file->guisection = _("Weather");

    opt.weather_coefficient_file = G_define_standard_option(G_OPT_F_INPUT);
    opt.weather_coefficient_file->key = "weather_coefficient_file";
    opt.weather_coefficient_file->label =
        _("Weather coefficients or weather coefficients mean values");
    opt.weather_coefficient_file->description =
        _("Input file with one weather coefficient map name per line");
    opt.weather_coefficient_file->required = NO;
    opt.weather_coefficient_file->guisection = _("Weather");

    opt.weather_coefficient_stddev_file = G_define_standard_option(G_OPT_F_INPUT);
    opt.weather_coefficient_stddev_file->key = "weather_coefficient_stddev_file";
    opt.weather_coefficient_stddev_file->label =
        _("Standard deviation of weather coefficient");
    opt.weather_coefficient_stddev_file->description =
        _("Input file with one standard deviation map name per line");
    opt.weather_coefficient_stddev_file->required = NO;
    opt.weather_coefficient_stddev_file->guisection = _("Weather");

    // Survival rate

    opt.survival_rate_file = G_define_standard_option(G_OPT_F_INPUT);
    opt.survival_rate_file->key = "survival_rate";
    opt.survival_rate_file->label =
        _("Input file with one survival rate raster map name per line");
    opt.survival_rate_file->description = _("Suvival rate is percentage (0-1)");
    opt.survival_rate_file->required = NO;
    opt.survival_rate_file->guisection = _("Weather");

    opt.survival_rate_month = G_define_option();
    opt.survival_rate_month->type = TYPE_INTEGER;
    opt.survival_rate_month->key = "survival_month";
    opt.survival_rate_month->label =
        _("Month when the pest or pathogen dies with given rate");
    opt.survival_rate_month->description =
        _("The survival rate is applied at the selected month and day");
    opt.survival_rate_month->required = NO;
    opt.survival_rate_month->guisection = _("Weather");

    opt.survival_rate_day = G_define_option();
    opt.survival_rate_day->type = TYPE_INTEGER;
    opt.survival_rate_day->key = "survival_day";
    opt.survival_rate_day->label =
        _("Day of selected month when the pest or pathogen dies with given rate");
    opt.survival_rate_day->description =
        _("The survival rate is applied at the selected month and day");
    opt.survival_rate_day->required = NO;
    opt.survival_rate_day->guisection = _("Weather");

    // Temperature and lethal temperature

    opt.lethal_temperature = G_define_option();
    opt.lethal_temperature->type = TYPE_DOUBLE;
    opt.lethal_temperature->key = "lethal_temperature";
    opt.lethal_temperature->label = _("Temperature at which the pest or pathogen dies");
    opt.lethal_temperature->description =
        _("The temerature unit must be the same as for the"
          "temerature raster map (typically degrees of Celsius)");
    opt.lethal_temperature->required = NO;
    opt.lethal_temperature->multiple = NO;
    opt.lethal_temperature->guisection = _("Weather");

    opt.lethal_temperature_months = G_define_option();
    opt.lethal_temperature_months->type = TYPE_INTEGER;
    opt.lethal_temperature_months->key = "lethal_month";
    opt.lethal_temperature_months->label =
        _("Month when the pest or pathogen dies due to low temperature");
    opt.lethal_temperature_months->description =
        _("The temperature unit must be the same as for the"
          "temperature raster map (typically degrees of Celsius)");
    // TODO: implement this as multiple
    opt.lethal_temperature_months->required = NO;
    opt.lethal_temperature_months->guisection = _("Weather");

    // TODO: rename coefs in interface and improve their descs
    opt.temperature_file = G_define_standard_option(G_OPT_F_INPUT);
    opt.temperature_file->key = "temperature_file";
    opt.temperature_file->label =
        _("Input file with one temperature raster map name per line");
    opt.temperature_file->description = _(
        "The temperature should be in actual temperature units (typically degrees of Celsius)");
    opt.temperature_file->required = NO;
    opt.temperature_file->guisection = _("Weather");

    opt.start_date = G_define_option();
    opt.start_date->type = TYPE_STRING;
    opt.start_date->key = "start_date";
    opt.start_date->description =
        _("Start date of the simulation in YYYY-MM-DD format");
    opt.start_date->required = YES;
    opt.start_date->guisection = _("Time");

    opt.end_date = G_define_option();
    opt.end_date->type = TYPE_STRING;
    opt.end_date->key = "end_date";
    opt.end_date->description = _("End date of the simulation in YYYY-MM-DD format");
    opt.end_date->required = YES;
    opt.end_date->guisection = _("Time");

    opt.seasonality = G_define_option();
    opt.seasonality->type = TYPE_STRING;
    opt.seasonality->key = "seasonality";
    opt.seasonality->label = _("Seasonal spread (from,to)");
    opt.seasonality->description =
        _("Spread limited to certain months (season), for example"
          " 5,9 for spread starting at the beginning of May and"
          " ending at the end of September");
    opt.seasonality->key_desc = "from,to";
    // opt.seasonality->options = "1-12";
    opt.seasonality->answer = const_cast<char*>("1,12");
    opt.seasonality->required = YES;
    opt.seasonality->multiple = NO;
    opt.seasonality->guisection = _("Time");

    opt.step_unit = G_define_option();
    opt.step_unit->type = TYPE_STRING;
    opt.step_unit->key = "step_unit";
    opt.step_unit->label = _("Unit of simulation steps");
    opt.step_unit->options = "day,week,month";
    opt.step_unit->answer = const_cast<char*>("month");
    opt.step_unit->descriptions = _(
        "day;Compute next simulation step every N days;week;Compute next simulation step every N weeks;month;Compute next simulation step every N months");
    opt.step_unit->required = YES;
    opt.step_unit->guisection = _("Time");

    opt.step_num_units = G_define_option();
    opt.step_num_units->type = TYPE_INTEGER;
    opt.step_num_units->key = "step_num_units";
    opt.step_num_units->answer = const_cast<char*>("1");
    opt.step_num_units->label = _("Number of days/weeks/months in each step");
    opt.step_num_units->description = _(
        "Step is given by number and unit, e.g. step_num_units=5 and step_unit=day means step is 5 days");
    opt.step_num_units->options = "1-";
    opt.step_num_units->required = YES;
    opt.step_num_units->guisection = _("Time");

    opt.output_frequency = G_define_option();
    opt.output_frequency->type = TYPE_STRING;
    opt.output_frequency->key = "output_frequency";
    opt.output_frequency->label = "Frequency of simulation output";
    opt.output_frequency->options = "yearly,monthly,weekly,daily,every_n_steps";
    opt.output_frequency->required = NO;
    opt.output_frequency->answer = const_cast<char*>("yearly");
    opt.output_frequency->guisection = _("Time");

    opt.output_frequency_n = G_define_option();
    opt.output_frequency_n->type = TYPE_INTEGER;
    opt.output_frequency_n->key = "output_frequency_n";
    opt.output_frequency_n->description = "Output frequency every N steps";
    opt.output_frequency_n->options = "1-";
    opt.output_frequency_n->answer = const_cast<char*>("1");
    opt.output_frequency_n->required = NO;
    opt.output_frequency_n->guisection = _("Time");

    opt.reproductive_rate = G_define_option();
    opt.reproductive_rate->type = TYPE_DOUBLE;
    opt.reproductive_rate->key = "reproductive_rate";
    opt.reproductive_rate->label =
        _("Number of spores or pest units produced by a single host");
    opt.reproductive_rate->description = _(
        "Number of spores or pest units produced by a single host under optimal weather conditions");
    opt.reproductive_rate->answer = const_cast<char*>("4.4");
    opt.reproductive_rate->guisection = _("Dispersal");

    opt.natural_kernel = G_define_option();
    opt.natural_kernel->type = TYPE_STRING;
    opt.natural_kernel->key = "natural_dispersal_kernel";
    opt.natural_kernel->label = _("Natural dispersal kernel type");
    opt.natural_kernel->answer = const_cast<char*>("cauchy");
    opt.natural_kernel->options = "cauchy,exponential";
    opt.natural_kernel->required = YES;
    opt.natural_kernel->guisection = _("Dispersal");

    opt.natural_scale = G_define_option();
    opt.natural_scale->type = TYPE_DOUBLE;
    opt.natural_scale->key = "natural_distance";
    opt.natural_scale->label = _("Distance parameter for natural dispersal kernel");
    opt.natural_scale->required = YES;
    opt.natural_scale->guisection = _("Dispersal");

    opt.natural_direction = G_define_option();
    opt.natural_direction->type = TYPE_STRING;
    opt.natural_direction->key = "natural_direction";
    opt.natural_direction->label = _("Direction of natural dispersal kernel");
    opt.natural_direction->description =
        _("Typically prevailing wind direction;"
          " none means that there is no directionality or no wind");
    opt.natural_direction->options = "N,NE,E,SE,S,SW,W,NW,NONE,none";
    opt.natural_direction->required = YES;
    opt.natural_direction->answer = const_cast<char*>("none");
    opt.natural_direction->guisection = _("Dispersal");

    opt.natural_kappa = G_define_option();
    opt.natural_kappa->type = TYPE_DOUBLE;
    opt.natural_kappa->key = "natural_direction_strength";
    opt.natural_kappa->label = _("Strength of direction of natural dispersal kernel");
    opt.natural_kappa->description =
        _("The kappa parameter of von Mises distribution"
          " (concentration);"
          " typically the strength of the wind direction");
    opt.natural_kappa->required = NO;
    opt.natural_kappa->guisection = _("Dispersal");

    opt.anthro_kernel = G_define_option();
    opt.anthro_kernel->type = TYPE_STRING;
    opt.anthro_kernel->key = "anthropogenic_dispersal_kernel";
    opt.anthro_kernel->label = _("Anthropogenic dispersal kernel type");
    opt.anthro_kernel->options = "cauchy,exponential";
    opt.anthro_kernel->guisection = _("Dispersal");

    opt.anthro_scale = G_define_option();
    opt.anthro_scale->type = TYPE_DOUBLE;
    opt.anthro_scale->key = "anthropogenic_distance";
    opt.anthro_scale->label =
        _("Distance parameter for anthropogenic dispersal kernel");
    opt.anthro_scale->guisection = _("Dispersal");

    opt.anthro_direction = G_define_option();
    opt.anthro_direction->type = TYPE_STRING;
    opt.anthro_direction->key = "anthropogenic_direction";
    opt.anthro_direction->label = _("Direction of anthropogenic dispersal kernel");
    opt.anthro_direction->description =
        _("Value none means that there is no directionality");
    opt.anthro_direction->options = "N,NE,E,SE,S,SW,W,NW,NONE,none";
    opt.anthro_direction->required = NO;
    opt.anthro_direction->answer = const_cast<char*>("none");
    opt.anthro_direction->guisection = _("Dispersal");

    opt.anthro_kappa = G_define_option();
    opt.anthro_kappa->type = TYPE_DOUBLE;
    opt.anthro_kappa->key = "anthropogenic_direction_strength";
    opt.anthro_kappa->label =
        _("Strength of direction of anthropogenic dispersal kernel");
    opt.anthro_kappa->description =
        _("The kappa parameter of von Mises distribution"
          " (concentration);"
          " typically the strength of the wind direction");
    opt.anthro_kappa->guisection = _("Dispersal");

    opt.percent_natural_dispersal = G_define_option();
    opt.percent_natural_dispersal->type = TYPE_DOUBLE;
    opt.percent_natural_dispersal->key = "percent_natural_dispersal";
    opt.percent_natural_dispersal->label = _("Percentage of natural dispersal");
    opt.percent_natural_dispersal->description =
        _("How often is the natural dispersal kernel used versus"
          " the anthropogenic dispersal kernel");
    opt.percent_natural_dispersal->options = "0-1";
    opt.percent_natural_dispersal->guisection = _("Dispersal");

    opt.dispersers_output = G_define_standard_option(G_OPT_R_OUTPUT);
    opt.dispersers_output->key = "dispersers_output";
    opt.dispersers_output->label = _("Output raster of disperses");
    opt.dispersers_output->description =
        _("Dispersers are accumulated over all steps and stochastic runs");
    opt.dispersers_output->required = NO;
    opt.dispersers_output->guisection = _("Output");

    opt.established_dispersers_output = G_define_standard_option(G_OPT_R_OUTPUT);
    opt.established_dispersers_output->key = "established_dispersers_output";
    opt.established_dispersers_output->label =
        _("Output raster of established disperses");
    opt.established_dispersers_output->description =
        _("Dispersers are accumulated over all steps and stochastic runs");
    opt.established_dispersers_output->required = NO;
    opt.established_dispersers_output->guisection = _("Output");

    opt.infected_to_dead_rate = G_define_option();
    opt.infected_to_dead_rate->type = TYPE_DOUBLE;
    opt.infected_to_dead_rate->key = "mortality_rate";
    opt.infected_to_dead_rate->label = _("Mortality rate of infected hosts");
    opt.infected_to_dead_rate->description =
        _("Percentage of infected hosts that die in a given year"
          " (hosts are removed from the infected pool)");
    opt.infected_to_dead_rate->options = "0-1";
    opt.infected_to_dead_rate->guisection = _("Mortality");

    opt.first_year_to_die = G_define_option();
    opt.first_year_to_die->type = TYPE_INTEGER;
    opt.first_year_to_die->key = "mortality_time_lag";
    opt.first_year_to_die->label =
        _("Time lag from infection until mortality can occur in years");
    opt.first_year_to_die->description =
        _("How many years it takes for an infected host to start dying"
          " (value 0 for hosts dying at the end of the first year)");
    opt.first_year_to_die->guisection = _("Mortality");

    opt.dead_series = G_define_standard_option(G_OPT_R_BASENAME_OUTPUT);
    opt.dead_series->key = "mortality_series";
    opt.dead_series->label = _("Basename for series of number of dead hosts");
    opt.dead_series->description =
        _("Basename for output series of number of dead hosts"
          " (requires mortality to be activated)");
    opt.dead_series->required = NO;
    opt.dead_series->guisection = _("Mortality");

    opt.mortality_frequency = G_define_option();
    opt.mortality_frequency->type = TYPE_STRING;
    opt.mortality_frequency->key = "mortality_frequency";
    opt.mortality_frequency->description = _("Mortality frequency");
    opt.mortality_frequency->options = "yearly,monthly,weekly,daily,every_n_steps";
    opt.mortality_frequency->required = NO;
    opt.mortality_frequency->guisection = _("Mortality");

    opt.mortality_frequency_n = G_define_option();
    opt.mortality_frequency_n->type = TYPE_INTEGER;
    opt.mortality_frequency_n->key = "mortality_frequency_n";
    opt.mortality_frequency_n->description = _("Mortality frequency every N steps");
    opt.mortality_frequency_n->options = "1-";
    opt.mortality_frequency_n->required = NO;
    opt.mortality_frequency_n->guisection = _("Mortality");

    flg.mortality = G_define_flag();
    flg.mortality->key = 'm';
    flg.mortality->label = _("Apply mortality");
    flg.mortality->description =
        _("After certain number of years, start removing dead hosts"
          " from the infected pool with a given rate");
    flg.mortality->guisection = _("Mortality");

    opt.seed = G_define_option();
    opt.seed->key = "random_seed";
    opt.seed->type = TYPE_INTEGER;
    opt.seed->required = NO;
    opt.seed->label = _("Seed for random number generator");
    opt.seed->description =
        _("The same seed can be used to obtain same results"
          " or random seed can be generated by other means.");
    opt.seed->guisection = _("Randomness");

    opt.seeds = G_define_option();
    opt.seeds->key = "random_seeds";
    opt.seeds->key_desc = "name=value";
    opt.seeds->type = TYPE_STRING;
    opt.seeds->required = NO;
    opt.seeds->multiple = YES;
    opt.seeds->label = _("Seeds for isolated random number generators");
    opt.seeds->description = _(
        "Multiple seeds for separate processes as a list of pairs key=value (comma-separated). "
        "Seeds must be provided for: disperser_generation,natural_dispersal,"
        "anthropogenic_dispersal,establishment,weather,movement,"
        "overpopulation,survival_rate,soil");
    opt.seeds->guisection = _("Randomness");

    flg.generate_seed = G_define_flag();
    flg.generate_seed->key = 's';
    flg.generate_seed->label = _("Generate random seed (result is non-deterministic)");
    flg.generate_seed->description =
        _("Automatically generates random seed for random number"
          " generator (use when you don't want to provide the seed option)");
    flg.generate_seed->guisection = _("Randomness");

    opt.runs = G_define_option();
    opt.runs->key = "runs";
    opt.runs->type = TYPE_INTEGER;
    opt.runs->required = NO;
    opt.runs->label = _("Number of simulation runs");
    opt.runs->description =
        _("The individual runs will obtain different seeds"
          " and will be averaged for the output");
    opt.runs->options = "1-";
    opt.runs->guisection = _("Randomness");

    opt.threads = G_define_option();
    opt.threads->key = "nprocs";
    opt.threads->type = TYPE_INTEGER;
    opt.threads->required = NO;
    opt.threads->description = _("Number of threads for parallel computing");
    opt.threads->options = "1-";
    opt.threads->guisection = _("Randomness");

    G_option_required(
        opt.average,
        opt.average_series,
        opt.single_series,
        opt.probability,
        opt.probability_series,
        opt.outside_spores,
        opt.stddev,
        opt.stddev_series,
        NULL);
    G_option_requires_all(opt.average_series, opt.output_frequency, NULL);
    G_option_requires_all(opt.single_series, opt.output_frequency, NULL);
    G_option_requires_all(opt.probability_series, opt.output_frequency, NULL);
    G_option_requires_all(opt.stddev_series, opt.output_frequency, NULL);
    G_option_exclusive(opt.seed, opt.seeds, flg.generate_seed, NULL);
    G_option_required(opt.seed, opt.seeds, flg.generate_seed, NULL);

    // weather
    G_option_collective(
        opt.moisture_coefficient_file, opt.temperature_coefficient_file, NULL);
    G_option_exclusive(
        opt.moisture_coefficient_file, opt.weather_coefficient_file, NULL);
    G_option_exclusive(
        opt.temperature_coefficient_file, opt.weather_coefficient_file, NULL);
    G_option_requires_all(
        opt.weather_coefficient_stddev_file, opt.weather_coefficient_file, NULL);

    // mortality
    // With flag, the lag, rate, and frequency are always required.
    // For simplicity of the code, mortality outputs are allowed only with the
    // corresponding main outputs for single run.
    G_option_collective(
        flg.mortality,
        opt.first_year_to_die,
        opt.infected_to_dead_rate,
        opt.mortality_frequency,
        NULL);
    G_option_requires_all(opt.dead_series, flg.mortality, opt.single_series, NULL);
    G_option_requires_all(opt.mortality_frequency_n, opt.mortality_frequency, NULL);
    // TODO: requires_all does not understand the default?
    // treatment_app needs to be removed from here and check separately
    G_option_requires_all(
        opt.treatments,
        opt.treatment_length,
        opt.treatment_date,
        opt.treatment_app,
        NULL);
    // lethal temperature options
    G_option_collective(
        opt.lethal_temperature,
        opt.lethal_temperature_months,
        opt.temperature_file,
        NULL);
    G_option_collective(opt.soil_survival_steps, opt.dispersers_to_soils, NULL);
    G_option_collective(
        opt.survival_rate_file, opt.survival_rate_month, opt.survival_rate_day, NULL);
    G_option_collective(opt.quarantine, opt.quarantine_output, NULL);

    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

    unsigned num_runs = 1;
    if (opt.runs->answer)
        num_runs = std::stoul(opt.runs->answer);

    unsigned threads = 1;
    if (opt.threads->answer)
        threads = std::stoul(opt.threads->answer);

    // check for file existence
    file_exists_or_fatal_error(opt.moisture_coefficient_file);
    file_exists_or_fatal_error(opt.temperature_coefficient_file);
    file_exists_or_fatal_error(opt.weather_coefficient_file);
    file_exists_or_fatal_error(opt.weather_coefficient_stddev_file);

    // Start creating the configuration.
    Config config;

    // model type
    config.model_type = opt.model_type->answer;
    ModelType model_type = model_type_from_string(config.model_type);
    if (model_type == ModelType::SusceptibleExposedInfected
        && !opt.latency_period->answer) {
        G_fatal_error(
            _("The option %s is required for %s=%s"),
            opt.latency_period->key,
            opt.model_type->key,
            opt.model_type->answer);
    }

    // get current computational region (for rows, cols and resolution)
    struct Cell_head window;
    G_get_window(&window);
    config.rows = window.rows;
    config.cols = window.cols;
    config.ew_res = window.ew_res;
    config.ns_res = window.ns_res;

    // Seasonality: Do you want the spread to be limited to certain months?
    if (!opt.seasonality->answer || opt.seasonality->answer[0] == '\0')
        G_fatal_error(_("The option %s cannot be empty"), opt.seasonality->key);
    seasonality_from_option(config, opt.seasonality);

    // set the spore rate
    config.reproductive_rate = std::stod(opt.reproductive_rate->answer);

    // TODO: how to support DispersalKernelType::None for natural_kernel?
    // perhaps the long-short should take the type instead of bool,
    // then T is anything else than None
    // TODO: should all kernels support None?
    config.natural_kernel_type = opt.natural_kernel->answer;
    config.natural_scale = std::stod(opt.natural_scale->answer);
    config.natural_direction = opt.natural_direction->answer;
    config.natural_kappa = 0;
    if (direction_from_string(config.natural_direction) != Direction::None
        && !opt.natural_kappa->answer)
        fatal_option_required_for_value(opt.natural_kappa, opt.natural_direction);
    else if (opt.natural_kappa->answer)
        config.natural_kappa = std::stod(opt.natural_kappa->answer);

    config.use_anthropogenic_kernel = false;
    if (opt.anthro_kernel->answer)
        config.anthro_kernel_type = opt.anthro_kernel->answer;
    if (kernel_type_from_string(config.anthro_kernel_type) != DispersalKernelType::None)
        config.use_anthropogenic_kernel = true;

    config.anthro_scale = 1;
    if (config.use_anthropogenic_kernel && !opt.anthro_scale->answer)
        fatal_option_required_for_value(opt.anthro_scale, opt.anthro_kernel);
    else if (opt.anthro_scale->answer)
        config.anthro_scale = std::stod(opt.anthro_scale->answer);

    // we allow none and an empty string
    config.anthro_direction = opt.anthro_direction->answer;

    config.anthro_kappa = 0;
    if (config.use_anthropogenic_kernel && !opt.anthro_kappa->answer)
        fatal_option_required_for_value(opt.anthro_kappa, opt.anthro_kernel);
    else if (opt.anthro_kappa->answer)
        config.anthro_kappa = std::stod(opt.anthro_kappa->answer);

    config.percent_natural_dispersal = 0.0;
    if (config.use_anthropogenic_kernel && !opt.percent_natural_dispersal->answer)
        fatal_option_required_for_value(
            opt.percent_natural_dispersal, opt.anthro_kernel);
    else if (opt.percent_natural_dispersal->answer)
        config.percent_natural_dispersal =
            std::stod(opt.percent_natural_dispersal->answer);

    // warn about limits to backwards compatibility
    // "none" is consistent with other GRASS GIS modules
    warn_about_depreciated_option_value(opt.natural_direction, "NONE", "none");
    warn_about_depreciated_option_value(opt.anthro_kernel, "NONE", "none");
    warn_about_depreciated_option_value(opt.anthro_direction, "NONE", "none");

    config.set_date_start(opt.start_date->answer);
    config.set_date_end(opt.end_date->answer);
    if (config.date_start() > config.date_end()) {
        G_fatal_error(_("Start date must precede the end date"));
    }

    config.set_step_unit(opt.step_unit->answer);
    config.set_step_num_units(std::stoi(opt.step_num_units->answer));

    config.output_frequency =
        opt.output_frequency->answer ? opt.output_frequency->answer : "";
    config.output_frequency_n =
        opt.output_frequency_n->answer ? std::stoi(opt.output_frequency_n->answer) : 0;

    config.latency_period_steps = 0;
    if (opt.latency_period->answer)
        config.latency_period_steps = std::stoi(opt.latency_period->answer);

    // mortality
    if (flg.mortality->answer) {
        config.use_mortality = true;
        config.mortality_time_lag =
            std::stoi(opt.first_year_to_die->answer);  // starts at 1 (same as the opt)
        config.mortality_rate = std::stod(opt.infected_to_dead_rate->answer);
        config.mortality_frequency = opt.mortality_frequency->answer;
        config.mortality_frequency_n =
            opt.mortality_frequency_n->answer
                ? std::stoi(opt.mortality_frequency_n->answer)
                : 0;
    }
    config.create_pest_host_table_from_parameters(1);
    std::vector<std::vector<double>> competency_table_data;
    competency_table_data.push_back({1, 1});
    competency_table_data.push_back({0, 0});
    config.read_competency_table(competency_table_data);

    if (opt.survival_rate_month->answer)
        config.survival_rate_month = std::stoi(opt.survival_rate_month->answer);
    if (opt.survival_rate_day->answer)
        config.survival_rate_day = std::stod(opt.survival_rate_day->answer);
    if (opt.survival_rate_file->answer)
        config.use_survival_rate = true;

    if (opt.lethal_temperature->answer)
        config.lethal_temperature = std::stod(opt.lethal_temperature->answer);
    if (opt.lethal_temperature_months->answer)
        config.lethal_temperature_month =
            std::stoi(opt.lethal_temperature_months->answer);
    if (opt.temperature_file->answer)
        config.use_lethal_temperature = true;

    config.use_spreadrates = false;
    if (opt.spread_rate_output->answer) {
        config.use_spreadrates = true;
        config.spreadrate_frequency = "yearly";
        config.spreadrate_frequency_n = 1;
    }
    config.use_quarantine = false;
    if (opt.quarantine_output->answer) {
        config.use_quarantine = true;
        config.quarantine_frequency = "yearly";
        config.quarantine_frequency_n = 1;
        if (opt.quarantine_directions->answer)
            config.quarantine_directions = opt.quarantine_directions->answer;
    }

    std::vector<string> moisture_names;
    std::vector<string> temperature_names;
    std::vector<string> weather_names;
    std::vector<string> weather_stddev_names;
    bool weather = false;
    bool moisture_temperature = false;
    bool weather_coefficient_distribution = false;
    if (opt.moisture_coefficient_file->answer
        && opt.temperature_coefficient_file->answer) {
        moisture_names = read_names(opt.moisture_coefficient_file->answer);
        temperature_names = read_names(opt.temperature_coefficient_file->answer);
        moisture_temperature = true;
    }
    if (opt.weather_coefficient_file->answer) {
        weather_names = read_names(opt.weather_coefficient_file->answer);
        weather = true;
    }
    if (opt.weather_coefficient_stddev_file->answer) {
        weather_stddev_names = read_names(opt.weather_coefficient_stddev_file->answer);
        weather_coefficient_distribution = true;
        if (weather_names.size() != weather_stddev_names.size()) {
            G_fatal_error(
                _("%s and %s must have the same size (not %d and %d)"),
                opt.weather_coefficient_file->key,
                opt.weather_coefficient_stddev_file->key,
                int(weather_names.size()),
                int(weather_stddev_names.size()));
        }
    }
    config.weather_size = weather_names.size();
    // Model gets pre-computed weather coefficient, so it does not
    // distinguish between these two.
    config.weather = weather || moisture_temperature;
    if (config.weather && weather_coefficient_distribution)
        config.weather_type = "probabilistic";
    else if (config.weather)
        config.weather_type = "deterministic";
    WeatherType weather_type = weather_type_from_string(config.weather_type);

    config.create_schedules();

    int num_mortality_steps = config.num_mortality_steps();
    if (flg.mortality->answer) {
        if (config.mortality_time_lag > num_mortality_steps) {
            G_fatal_error(
                _("%s is too large (%d). It must be smaller or "
                  " equal than number of simulation years (%d)."),
                opt.first_year_to_die->key,
                config.mortality_time_lag,
                num_mortality_steps);
        }
    }
    else
        num_mortality_steps = 1;

    if (opt.seed->answer) {
        config.random_seed = std::stoul(opt.seed->answer);
        G_verbose_message(
            _("Using random seed from %s option: %u"),
            opt.seed->key,
            config.random_seed);
    }
    else if (opt.seeds->answer) {
        config.read_seeds(opt.seeds->answer, ',', '=');
        try {
            validate_random_number_generator_provider_config(config);
        }
        catch (const std::invalid_argument& error) {
            G_fatal_error(
                _("%s is incomplete or incorrectly formatted: %s"),
                opt.seeds->key,
                error.what());
        }

        G_verbose_message(
            _("Using random seeds from %s option: %s"),
            opt.seeds->key,
            opt.seeds->answer);
    }
    else {
        // Flag or option is required, so no further check is needed here.
        // getting random seed using GRASS library
        // std::random_device is deterministic in MinGW (#338)
        config.random_seed = G_srand48_auto();
        G_verbose_message(
            _("Generated random seed (-%c): %u"),
            flg.generate_seed->key,
            config.random_seed);
    }

    // read the suspectible UMCA raster image
    Img species_rast = raster_from_grass_integer(opt.host->answer);

    // read the living trees raster image
    Img lvtree_rast = raster_from_grass_integer(opt.total_plants->answer);

    // read the initial infected oaks image
    Img I_species_rast = raster_from_grass_integer(opt.infected->answer);

    // create the initial suspectible oaks image
    Img S_species_rast = species_rast - I_species_rast;

    // Indices of suitable cells (i, j of all hosts)
    std::vector<std::vector<int>> suitable_cells = get_suitable_cells(species_rast);

    std::vector<DImg> survival_rates;
    if (opt.survival_rate_file->answer) {
        survival_rates =
            file_to_series(opt.survival_rate_file, config.num_survival_rate());
    }

    std::vector<DImg> actual_temperatures;
    if (opt.temperature_file->answer) {
        actual_temperatures = file_to_series(opt.temperature_file, config.num_lethal());
    }

    std::vector<DImg> weather_coefficients;
    if (weather || moisture_temperature)
        weather_coefficients.resize(config.scheduler().get_num_steps());
    std::vector<DImg> weather_coefficient_stddevs;
    if (weather_coefficient_distribution)
        weather_coefficient_stddevs.resize(config.scheduler().get_num_steps());

    bool use_soils = false;
    int soil_survival_steps = 0;
    if (opt.soil_survival_steps->answer) {
        use_soils = true;
        soil_survival_steps = std::stoi(opt.soil_survival_steps->answer);
        config.dispersers_to_soils_percentage =
            std::stod(opt.dispersers_to_soils->answer);
    }

    using SpreadModel = Model<Img, DImg, DImg::IndexType>;

    // treatments
    if (get_num_answers(opt.treatments) != get_num_answers(opt.treatment_date)
        && get_num_answers(opt.treatment_date)
               != get_num_answers(opt.treatment_length)) {
        G_fatal_error(
            _("%s=, %s= and %s= must have the same number of values"),
            opt.treatments->key,
            opt.treatment_date->key,
            opt.treatment_length->key);
    }
    // the default here should be never used
    TreatmentApplication treatment_app = TreatmentApplication::Ratio;
    if (opt.treatment_app->answer)
        treatment_app = treatment_app_enum_from_string(opt.treatment_app->answer);
    Treatments<SpreadModel::StandardSingleHostPool, DImg> treatments(
        config.scheduler());
    config.use_treatments = false;
    if (opt.treatments->answers) {
        for (int i_t = 0; opt.treatment_date->answers[i_t]; i_t++) {
            DImg tr = raster_from_grass_float(opt.treatments->answers[i_t]);
            treatments.add_treatment(
                tr,
                treatment_date_from_string(opt.treatment_date->answers[i_t]),
                std::stoul(opt.treatment_length->answers[i_t]),
                treatment_app);
            config.use_treatments = true;
        }
    }

    // build the model object
    std::vector<SpreadModel> models;
    std::vector<Img> dispersers;
    std::vector<Img> established_dispersers;
    std::vector<Img> sus_species_rasts(num_runs, S_species_rast);
    std::vector<Img> inf_species_rasts(num_runs, I_species_rast);
    std::vector<Img> total_species_rasts(num_runs, species_rast);
    std::vector<Img> resistant_rasts(num_runs, Img(S_species_rast, 0));
    std::vector<Img> total_exposed_rasts(num_runs, Img(S_species_rast, 0));
    std::vector<SpreadModel::StandardMultiHostPool> multi_host_pools;
    std::vector<PestPool<Img, DImg, int>> pest_pools;

    // We always create at least one exposed for simplicity, but we
    // could also just leave it empty.
    std::vector<std::vector<Img>> exposed_vectors(
        num_runs,
        std::vector<Img>(
            config.latency_period_steps + 1,
            Img(S_species_rast.rows(), S_species_rast.cols(), 0)));

    // infected cohort for each year (index is cohort age)
    // age starts with 0 (in year 1), 0 is oldest
    std::vector<std::vector<Img>> mortality_tracker_vector(
        num_runs, std::vector<Img>(num_mortality_steps, Img(S_species_rast, 0)));

    // we are using only the first dead img for visualization, but for
    // parallelization we need all allocated anyway
    std::vector<Img> dead_in_current_year(num_runs, Img(S_species_rast, 0));
    // dead trees accumulated over years
    // TODO: allow only when series as single run
    Img accumulated_dead(Img(S_species_rast, 0));

    models.reserve(num_runs);
    multi_host_pools.reserve(num_runs);
    pest_pools.reserve(num_runs);
    dispersers.reserve(num_runs);
    established_dispersers.reserve(num_runs);
    auto tmp_seeds = config.random_seeds;
    auto tmp_seed_value = config.random_seed;
    for (unsigned i = 0; i < num_runs; ++i) {
        Config config_copy = config;
        if (opt.seeds->answer) {
            for (auto& item : tmp_seeds) {
                config_copy.random_seeds[item.first] = item.second++;
            }
        }
        else {
            config_copy.random_seed = tmp_seed_value++;
        }
        try {
            models.emplace_back(config_copy);
        }
        catch (const std::invalid_argument& error) {
            G_fatal_error(_("Model configuration is invalid: %s"), error.what());
        }
        // For the model itself, this does not need to be initialized to 0 because the
        // model only sets and then looks at suitable cells. However, the output is
        // taken from all cells, so even the untouched cells need a value.
        dispersers.emplace_back(I_species_rast.rows(), I_species_rast.cols(), 0);
        established_dispersers.emplace_back(
            I_species_rast.rows(), I_species_rast.cols(), 0);
    }
    std::vector<Img> dispersers_rasts(num_runs, Img(S_species_rast, 0));
    std::vector<Img> established_dispersers_rasts(num_runs, Img(S_species_rast, 0));

    std::vector<std::vector<Img>> soil_reservoirs(
        use_soils ? num_runs : 0,
        std::vector<Img>(
            use_soils ? soil_survival_steps : 0,
            Img(I_species_rast.rows(), I_species_rast.cols(), 0)));
    if (use_soils) {
        for (unsigned i = 0; i < num_runs; ++i) {
            models[i].activate_soils(soil_reservoirs[i]);
        }
    }
    std::vector<std::vector<std::tuple<int, int>>> outside_spores(num_runs);

    // One host pool for each run (not multi-host case).
    std::vector<std::unique_ptr<SpreadModel::StandardSingleHostPool>> host_pools;
    host_pools.reserve(num_runs);
    std::vector<
        std::unique_ptr<pops::PestHostTable<SpreadModel::StandardSingleHostPool>>>
        pest_host_tables;
    pest_host_tables.reserve(num_runs);
    std::vector<
        std::unique_ptr<pops::CompetencyTable<SpreadModel::StandardSingleHostPool>>>
        competency_tables;
    competency_tables.reserve(num_runs);

    for (unsigned run = 0; run < num_runs; ++run) {
        pest_host_tables.emplace_back(
            new pops::PestHostTable<SpreadModel::StandardSingleHostPool>(
                config, models[run].environment()));
        competency_tables.emplace_back(
            new pops::CompetencyTable<SpreadModel::StandardSingleHostPool>(
                config, models[run].environment()));

        host_pools.emplace_back(new SpreadModel::StandardSingleHostPool(
            model_type_from_string(config.model_type),
            sus_species_rasts[run],
            exposed_vectors[run],
            config.latency_period_steps,
            inf_species_rasts[run],
            total_exposed_rasts[run],
            resistant_rasts[run],
            mortality_tracker_vector[run],
            dead_in_current_year[run],
            total_species_rasts[run],
            models[run].environment(),
            config.generate_stochasticity,
            config.reproductive_rate,
            config.establishment_stochasticity,
            config.establishment_probability,
            config.rows,
            config.cols,
            suitable_cells));
        std::vector<SpreadModel::StandardSingleHostPool*> tmp = {host_pools[run].get()};
        multi_host_pools.emplace_back(tmp, config);
        multi_host_pools[run].set_pest_host_table(*pest_host_tables[run]);
        multi_host_pools[run].set_competency_table(*competency_tables[run]);

        pest_pools.emplace_back(
            dispersers[run], established_dispersers[run], outside_spores[run]);
    }
    std::vector<SpreadRateAction<SpreadModel::StandardMultiHostPool, int>> spread_rates(
        num_runs,
        SpreadRateAction<SpreadModel::StandardMultiHostPool, int>(
            multi_host_pools[0],
            config.rows,
            config.cols,
            config.ew_res,
            config.ns_res,
            config.use_spreadrates ? config.rate_num_steps() : 0));
    // Quarantine escape tracking
    Img quarantine_rast(S_species_rast, 0);
    if (config.use_quarantine)
        quarantine_rast = raster_from_grass_integer(opt.quarantine->answer);
    std::vector<QuarantineEscapeAction<Img>> escape_infos(
        num_runs,
        QuarantineEscapeAction<Img>(
            quarantine_rast,
            window.ew_res,
            window.ns_res,
            config.use_quarantine ? config.quarantine_num_steps() : 0,
            config.quarantine_directions));

    // Unused movements
    std::vector<std::vector<int>> movements;

    std::vector<unsigned> unresolved_steps;
    unresolved_steps.reserve(config.scheduler().get_num_steps());

    G_verbose_message(_("Number of steps: %u"), config.scheduler().get_num_steps());

    // main simulation loop
    unsigned current_index = 0;
    for (; current_index < config.scheduler().get_num_steps(); ++current_index) {
        unresolved_steps.push_back(current_index);

        // if all the hosts are infected, then exit
        if (all_infected(S_species_rast)) {
            G_warning(
                "In step %d all suspectible hosts are infected, ending simulation.",
                current_index);
            break;
        }

        // check whether the spore occurs in the month
        // At the end of the year, run simulation for all unresolved
        // steps in one chunk.
        if (config.output_schedule()[current_index]
            || current_index == config.scheduler().get_num_steps() - 1) {
            unsigned step_in_chunk = 0;
            // get weather for all the steps in chunk
            if (weather) {
                for (auto step : unresolved_steps) {
                    auto weather_step = config.simulation_step_to_weather_step(step);
                    if (moisture_temperature) {
                        DImg moisture(
                            raster_from_grass_float(moisture_names[weather_step]));
                        DImg temperature(
                            raster_from_grass_float(temperature_names[weather_step]));
                        weather_coefficients[step_in_chunk] = moisture * temperature;
                    }
                    else if (weather_type == WeatherType::Probabilistic) {
                        weather_coefficients[step_in_chunk] =
                            raster_from_grass_float(weather_names[weather_step]);
                        weather_coefficient_stddevs[step_in_chunk] =
                            raster_from_grass_float(weather_stddev_names[weather_step]);
                    }
                    else {
                        weather_coefficients[step_in_chunk] =
                            raster_from_grass_float(weather_names[weather_step]);
                    }
                    ++step_in_chunk;
                }
            }

// stochastic simulation runs
#pragma omp parallel for num_threads(threads)
            for (unsigned run = 0; run < num_runs; run++) {
                // actual runs of the simulation for each step
                int weather_step = 0;
                for (auto step : unresolved_steps) {
                    dead_in_current_year[run].zero();
                    if (weather_type == WeatherType::Probabilistic) {
                        models[run].environment().update_weather_from_distribution(
                            weather_coefficients[weather_step],
                            weather_coefficient_stddevs[weather_step],
                            models[run].random_number_generator());
                    }
                    else if (weather_type == WeatherType::Deterministic) {
                        models[run].environment().update_weather_coefficient(
                            weather_coefficients[weather_step]);
                    }
                    models[run].run_step(
                        step,
                        multi_host_pools[run],
                        pest_pools[run],
                        lvtree_rast,
                        treatments,
                        actual_temperatures,
                        survival_rates,
                        spread_rates[run],
                        escape_infos[run],
                        quarantine_rast,
                        movements,
                        Network<Img::IndexType>::null_network());
                    ++weather_step;
                    if (opt.dispersers_output->answer)
                        dispersers_rasts[run] += dispersers[run];
                    if (opt.established_dispersers_output->answer)
                        established_dispersers_rasts[run] +=
                            established_dispersers[run];
                }
            }

            unresolved_steps.clear();
            if (config.output_schedule()[current_index]) {
                // output
                Step interval = config.scheduler().get_step(current_index);
                if (opt.single_series->answer) {
                    string name =
                        generate_name(opt.single_series->answer, interval.end_date());
                    raster_to_grass(
                        inf_species_rasts[0],
                        name,
                        "Occurrence from a single stochastic run",
                        interval.end_date());
                }
                if ((opt.average_series->answer) || opt.stddev_series->answer) {
                    // aggregate in the series
                    DImg average_raster(
                        I_species_rast.rows(), I_species_rast.cols(), 0);
                    average_raster.zero();
                    for (unsigned i = 0; i < num_runs; i++)
                        average_raster += inf_species_rasts[i];
                    average_raster /= num_runs;
                    if (opt.average_series->answer) {
                        // write result
                        // date is always end of the year, even for seasonal spread
                        string name = generate_name(
                            opt.average_series->answer, interval.end_date());
                        raster_to_grass(
                            average_raster,
                            name,
                            "Average occurrence from all stochastic runs",
                            interval.end_date());
                        write_average_area(
                            inf_species_rasts,
                            name.c_str(),
                            window.ew_res,
                            window.ns_res,
                            suitable_cells);
                    }
                    if (opt.stddev_series->answer) {
                        DImg stddev(I_species_rast.rows(), I_species_rast.cols(), 0);
                        for (unsigned i = 0; i < num_runs; i++) {
                            auto tmp = inf_species_rasts[i] - average_raster;
                            stddev += tmp * tmp;
                        }
                        stddev /= num_runs;
                        stddev.for_each([](Float& a) { a = std::sqrt(a); });
                        string name = generate_name(
                            opt.stddev_series->answer, interval.end_date());
                        string title =
                            "Standard deviation of average"
                            " occurrence from all stochastic runs";
                        raster_to_grass(stddev, name, title, interval.end_date());
                    }
                }
                if (opt.probability_series->answer) {
                    DImg probability(I_species_rast.rows(), I_species_rast.cols(), 0);
                    for (unsigned i = 0; i < num_runs; i++) {
                        Img tmp = inf_species_rasts[i];
                        tmp.for_each([](Integer& a) { a = bool(a); });
                        probability += tmp;
                    }
                    probability *= 100;  // prob from 0 to 100
                    probability /= num_runs;
                    string name = generate_name(
                        opt.probability_series->answer, interval.end_date());
                    string title = "Probability of occurrence";
                    raster_to_grass(probability, name, title, interval.end_date());
                }
                if (config.use_mortality && opt.dead_series->answer) {
                    accumulated_dead += dead_in_current_year[0];
                    if (opt.dead_series->answer) {
                        string name =
                            generate_name(opt.dead_series->answer, interval.end_date());
                        raster_to_grass(
                            accumulated_dead,
                            name,
                            "Number of dead hosts to date",
                            interval.end_date());
                    }
                }
            }
        }
    }
    Step interval = config.scheduler().get_step(--current_index);
    if (opt.average->answer || opt.stddev->answer) {
        // aggregate
        DImg average_raster(I_species_rast.rows(), I_species_rast.cols(), 0);
        for (unsigned i = 0; i < num_runs; i++)
            average_raster += inf_species_rasts[i];
        average_raster /= num_runs;
        if (opt.average->answer) {
            // write final result
            raster_to_grass(
                average_raster,
                opt.average->answer,
                "Average occurrence from all stochastic runs",
                interval.end_date());
            write_average_area(
                inf_species_rasts,
                opt.average->answer,
                window.ew_res,
                window.ns_res,
                suitable_cells);
        }
        if (opt.stddev->answer) {
            DImg stddev(average_raster.rows(), average_raster.cols(), 0);
            for (unsigned i = 0; i < num_runs; i++) {
                auto tmp = inf_species_rasts[i] - average_raster;
                stddev += tmp * tmp;
            }
            stddev /= num_runs;
            stddev.for_each([](Float& a) { a = std::sqrt(a); });
            raster_to_grass(
                stddev,
                opt.stddev->answer,
                opt.stddev->description,
                interval.end_date());
        }
    }
    if (opt.probability->answer) {
        DImg probability(I_species_rast.rows(), I_species_rast.cols(), 0);
        for (unsigned i = 0; i < num_runs; i++) {
            Img tmp = inf_species_rasts[i];
            tmp.for_each([](Integer& a) { a = bool(a); });
            probability += tmp;
        }
        probability *= 100;  // prob from 0 to 100
        probability /= num_runs;
        raster_to_grass(
            probability,
            opt.probability->answer,
            "Probability of occurrence",
            interval.end_date());
    }
    if (opt.outside_spores->answer) {
        Cell_head region;
        Rast_get_window(&region);
        struct Map_info Map;
        struct line_pnts* Points;
        struct line_cats* Cats;
        if (Vect_open_new(&Map, opt.outside_spores->answer, WITHOUT_Z) < 0)
            G_fatal_error(
                _("Unable to create vector map <%s>"), opt.outside_spores->answer);

        Points = Vect_new_line_struct();
        Cats = Vect_new_cats_struct();

        for (unsigned i = 0; i < num_runs; i++) {
            for (unsigned j = 0; j < outside_spores[i].size(); j++) {
                int row = std::get<0>(outside_spores[i][j]);
                int col = std::get<1>(outside_spores[i][j]);
                double n = Rast_row_to_northing(row, &region);
                double e = Rast_col_to_easting(col, &region);
                Vect_reset_line(Points);
                Vect_reset_cats(Cats);
                Vect_append_point(Points, e, n, 0);
                Vect_cat_set(Cats, 1, i + 1);
                Vect_write_line(&Map, GV_POINT, Points, Cats);
            }
        }
        Vect_hist_command(&Map);
        Vect_set_map_name(&Map, "Dispersers escaped outside computational region");
        Vect_write_header(&Map);
        Vect_build(&Map);
        Vect_close(&Map);
        struct TimeStamp timestamp;
        date_to_grass(interval.end_date(), &timestamp);
        G_write_vector_timestamp(opt.outside_spores->answer, NULL, &timestamp);
        Vect_destroy_line_struct(Points);
        Vect_destroy_cats_struct(Cats);
    }
    if (opt.spread_rate_output->answer) {
        FILE* fp = G_open_option_file(opt.spread_rate_output);
        fprintf(fp, "year,N,S,E,W\n");
        for (unsigned step = 0; step < config.scheduler().get_num_steps(); step++) {
            double n, s, e, w;
            if (config.spread_rate_schedule()[step]) {
                unsigned i =
                    simulation_step_to_action_step(config.spread_rate_schedule(), step);
                std::tie(n, s, e, w) = average_spread_rate(spread_rates, i);
                int year = config.scheduler().get_step(step).end_date().year();
                fprintf(
                    fp,
                    "%d,%.0f,%.0f,%.0f,%.0f\n",
                    year,
                    isnan(n) ? n : round(n),
                    isnan(s) ? s : round(s),
                    isnan(e) ? e : round(e),
                    isnan(w) ? w : round(w));
            }
        }
        G_close_option_file(fp);
    }
    if (opt.quarantine_output->answer) {
        FILE* fp = G_open_option_file(opt.quarantine_output);
        std::string output =
            write_quarantine_escape(escape_infos, config.quarantine_num_steps());
        fprintf(fp, "%s", output.c_str());
        G_close_option_file(fp);
    }
    if (opt.dispersers_output->answer) {
        Img dispersers_rasts_sum(I_species_rast.rows(), I_species_rast.cols(), 0);
        for (unsigned i = 0; i < num_runs; i++)
            dispersers_rasts_sum += dispersers_rasts[i];
        if (opt.dispersers_output->answer) {
            // write final result
            raster_to_grass(
                dispersers_rasts_sum,
                opt.dispersers_output->answer,
                "Sum of all dispersers",
                interval.end_date());
        }
    }
    if (opt.established_dispersers_output->answer) {
        Img established_dispersers_rasts_sum(
            I_species_rast.rows(), I_species_rast.cols(), 0);
        for (unsigned i = 0; i < num_runs; i++)
            established_dispersers_rasts_sum += established_dispersers_rasts[i];
        if (opt.established_dispersers_output->answer) {
            // write final result
            raster_to_grass(
                established_dispersers_rasts_sum,
                opt.established_dispersers_output->answer,
                "Sum of all established dispersers",
                interval.end_date());
        }
    }
    return 0;
}
