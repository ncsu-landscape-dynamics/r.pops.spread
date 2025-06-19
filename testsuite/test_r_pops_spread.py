#!/usr/bin/env python3

"""Test of r.pops.spread.

.. moduleauthor:: Anna Petrasova
.. moduleauthor:: Vaclav Petras
"""

import concurrent.futures
import tempfile
from pathlib import Path

import grass.script as gs
from grass.gunittest.case import TestCase
from grass.gunittest.main import test
from grass.gunittest.gmodules import call_module


def items_to_file(items, filename):
    """Save list of items to a file, one item per line"""
    with open(filename, mode="w", encoding="utf-8") as file:
        file.write("\n".join(items))


# Tests use dict function call because the syntax is simpler and there is
# a lot values provided like that.
# pylint: disable=use-dict-literal
# We have one setup and many tests.
# pylint: disable=too-many-lines,too-many-public-methods


class TestSpread(TestCase):
    """Tests of r.pops.spread"""

    @classmethod
    def setUpClass(cls):
        """Create input data from the full NC SPM dataset"""
        cls.use_temp_region()

        # Use simulation resolution.
        cls.runModule("g.region", raster="lsat7_2002_30", res=85.5, flags="a")
        cls.runModule(
            "r.mapcalc",
            expression=(
                "ndvi = double(lsat7_2002_40 - lsat7_2002_30)"
                " / double(lsat7_2002_40 + lsat7_2002_30)"
            ),
        )
        cls.runModule(
            "r.mapcalc",
            expression="host = round(if(ndvi > 0, graph(ndvi, 0, 0, 1, 20), 0))",
        )
        cls.runModule(
            "r.mapcalc", expression="host_nulls = if(host == 0, null(), host)"
        )
        cls.runModule(
            "v.to.rast", input="railroads", output="infection_", use="val", value=1
        )
        cls.runModule("r.null", map="infection_", null=0)
        cls.runModule("r.mapcalc", expression="infection = if(host > 0, infection_, 0)")
        cls.runModule(
            "r.mapcalc",
            expression="infection_nulls = if(infection == 0, null(), infection)",
        )
        cls.runModule("r.mapcalc", expression="max_host = 100")
        cls.runModule(
            "r.circle",
            flags="b",
            output="circle",
            coordinates=[639445, 218237],
            max=2000,
        )
        cls.runModule("r.mapcalc", expression="treatment = if(isnull(circle), 0, 0.8)")
        cls.runModule("r.mapcalc", expression="one = 1")
        cls.runModule("r.mapcalc", expression="zero = 0")

        # We can remove the directory only at the end, so we can't use with statement.
        cls.tmp_dir = (
            tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        )

        years_2019_2022 = 4
        months_in_year = 12
        weeks_in_year = 52

        cls.steps_in_2019_2022_monthly = years_2019_2022 * months_in_year
        # File with 100% survival rate time series
        cls.survival_rate_100_percent = str(
            Path(cls.tmp_dir.name) / "survival_rate_100_percent.txt"
        )
        items_to_file(
            cls.steps_in_2019_2022_monthly * ["one"], cls.survival_rate_100_percent
        )
        # File with 0% survival rate time series
        cls.survival_rate_0_percent = str(
            Path(cls.tmp_dir.name) / "survival_rate_0_percent.txt"
        )
        items_to_file(
            cls.steps_in_2019_2022_monthly * ["zero"], cls.survival_rate_0_percent
        )

        # Use weather resolution (extent will be larger than simulation extent).
        cls.runModule("g.region", raster="lsat7_2002_30", res=1000, flags="a")

        cls.weather_names = [f"weather_{i}" for i in range(weeks_in_year)]
        cls.weather_stddev_names = [f"weather_stddev_{i}" for i in range(weeks_in_year)]

        def generate_random_raster(name, low, high):
            gs.mapcalc(
                (f"{name} = rand(double({low}), double({high}))"),
                seed=1,
                superquiet=True,
            )
            return name

        # Number of processes based on the machine.
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(generate_random_raster, name, low=0, high=1)
                for name in cls.weather_names
            ]
            futures.extend(
                executor.submit(generate_random_raster, name, low=0, high=0)
                for name in cls.weather_stddev_names
            )
            concurrent.futures.wait(futures)
        cls.weather_file = str(Path(cls.tmp_dir.name) / "weather_coefficient.txt")
        items_to_file(years_2019_2022 * cls.weather_names, cls.weather_file)

        cls.weather_stddev_file = str(
            Path(cls.tmp_dir.name) / "weather_coefficient_stddev.txt"
        )
        items_to_file(
            (years_2019_2022) * cls.weather_stddev_names, cls.weather_stddev_file
        )

        # The zero raster was produced earlier.
        cls.weather_zero_stddev_file = str(
            Path(cls.tmp_dir.name) / "weather_coefficient_zero_stddev.txt"
        )
        items_to_file(
            (years_2019_2022 * weeks_in_year) * ["zero"], cls.weather_zero_stddev_file
        )

        # Quarantine
        cls.runModule("g.region", raster="lsat7_2002_30", res=85.5, flags="a")
        cls.runModule(
            "r.circle",
            output="raw_infected_patch",
            coordinates=[639300, 220900],
            max=100,
            flags="b",
        )
        gs.mapcalc("infected_patch = min(raw_infected_patch, host)", superquiet=True)
        cls.runModule(
            "g.region",
            n="n-800",
            s="s+800",
            e="e-800",
            w="w+800",
            align="lsat7_2002_30",
        )
        gs.mapcalc("quarantine = 1", superquiet=True)

        # Use simulation resolution for the computations.
        cls.runModule("g.region", raster="lsat7_2002_30", res=85.5, flags="a")

    @classmethod
    def tearDownClass(cls):
        """Remove the generated input data"""
        cls.tmp_dir.cleanup()
        cls.del_temp_region()
        cls.runModule(
            "g.remove",
            flags="f",
            type="raster",
            name=[
                "max_host",
                "infection_",
                "infection",
                "infection_nulls",
                "host",
                "host_nulls",
                "ndvi",
                "circle",
                "treatment",
                "one",
                "zero",
                "quarantine",
                "infected_patch",
                "raw_infected_patch",
                *cls.weather_names,
                *cls.weather_stddev_names,
            ],
        )

    def tearDown(self):
        """Remove maps after each test method"""
        self.runModule(
            "g.remove",
            flags="f",
            type=["raster", "vector"],
            pattern="average*,single*,stddev*,probability*,dead*,*dispersers",
        )

    def test_outputs(self):
        """Check standard outputs including last item in the time series"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")
        end = end[:4]
        self.assertRasterExists(f"average_{end}_12_31")
        self.assertRasterExists(f"probability_{end}_12_31")
        self.assertRasterExists(f"single_{end}_12_31")
        self.assertRasterExists(f"stddev_{end}_12_31")

        ref_float = dict(datatype="DCELL")
        ref_int = dict(datatype="CELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)
        self.assertRasterFitsInfo(raster=f"single_{end}_12_31", reference=ref_int)
        self.assertRasterFitsInfo(raster=f"average_{end}_12_31", reference=ref_float)
        self.assertRasterFitsInfo(
            raster=f"probability_{end}_12_31", reference=ref_float
        )
        self.assertRasterFitsInfo(raster=f"stddev_{end}_12_31", reference=ref_float)

        values = dict(null_cells=0, min=0, max=18, mean=1.777)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=33.664)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=7.547, mean=0.945)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

    def test_quarantine(self):
        """Check quarantine output"""
        start = "2019-01-01"
        end = "2022-12-31"
        output = call_module(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infected_patch",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            quarantine="quarantine",
            quarantine_output="-",
            output_frequency="yearly",
            random_seed=1,
            runs=5,
            nprocs=5,
        )
        reference = (
            "step,escape_probability,dist0,dir0,dist1,dir1,dist2,dir2,dist3,dir3,dist4,dir4\n"
            "0,0.0,4873.0,90,4018.0,180,4873.0,90,4873.0,90,4873.0,90\n"
            "1,0.0,4873.0,90,2052.0,180,4873.0,90,2394.0,180,4873.0,90\n"
            "2,0.4,1966.0,180,,,,,171.0,90,4873.0,90\n"
            "3,0.8,,,,,,,,,2137.0,90"
        )
        self.assertMultiLineEqual(output.strip(), reference)

    def test_weather_deterministic_strict(self):
        """Check deterministic weather with close-to-exact global statistics"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            weather_coefficient_file=self.weather_file,
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
        )
        test_date = "2021_12_31"
        end_year = end[:4]

        # Final outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.163)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=4.124)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=5.879, mean=0.114)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

        # Time-series outputs
        values = dict(null_cells=0, min=0, max=18.0, mean=0.098)
        self.assertRasterFitsUnivar(
            raster=f"average_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=100, mean=2.777)
        self.assertRasterFitsUnivar(
            raster=f"probability_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=5.643, mean=0.067)
        self.assertRasterFitsUnivar(
            raster=f"stddev_{test_date}", reference=values, precision=0.001
        )

        # Single run outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.110)
        self.assertRasterFitsUnivar(
            raster=f"single_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=18, mean=0.171)
        self.assertRasterFitsUnivar(
            raster=f"single_{end_year}_12_31", reference=values, precision=0.001
        )

    def test_weather_deterministic_lenient(self):
        """Check deterministic weather with global statistics and large tolerance"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            weather_coefficient_file=self.weather_file,
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
            overwrite=True,
        )
        test_date = "2021_12_31"
        end_year = end[:4]

        # The reference values were obtained from a run with 100 stochastic runs
        # with seed 1 and 2-3 non-zero digits were kept. The precision was chosen
        # so additional 300 runs with 1 stochastic run each and different seeds than
        # the original set (101-400) would each still pass the test.
        # This was done by modifying the test code to loop over a range of seeds with
        # random_seed=seed, and overwrite=True and increasing the required precision
        # (that is tolerance) as needed.
        # In case precision would become irrelevant for the smaller value (mean), the
        # test was split into two (the required precision is absolute, not relative).
        # Additional 10 runs with 100 stochastic runs each and different seeds
        # (401, 501, ..., 1301) also pass the test after small adjustments of required
        # precision. Additional changes for stddev were needed for seed 1 and number
        # of runs 5.

        # Final outputs
        values = dict(null_cells=0, min=0, max=18)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=1)
        values = dict(mean=0.170)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.05)
        values = dict(null_cells=0, min=0, max=100, mean=4.6)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=1.4
        )
        values = dict(null_cells=0, min=0, max=5.46)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.46)
        values = dict(mean=0.25)
        # For single run stddev is zero, but even for many runs the tolerance needs to be
        # almost 100% of the value.
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.25)

        # Time-series outputs
        values = dict(null_cells=0, min=0, max=17.9, mean=0.10)
        self.assertRasterFitsUnivar(
            raster=f"average_{test_date}", reference=values, precision=2.9
        )
        values = dict(null_cells=0, min=0, max=100, mean=3.03)
        self.assertRasterFitsUnivar(
            raster=f"probability_{test_date}", reference=values, precision=0.5
        )
        values = dict(null_cells=0, min=0, max=4.98)
        self.assertRasterFitsUnivar(
            raster=f"stddev_{test_date}", reference=values, precision=0.75
        )
        values = dict(mean=0.14)
        self.assertRasterFitsUnivar(
            raster=f"stddev_{test_date}", reference=values, precision=0.26
        )

        # Single run outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.10)
        self.assertRasterFitsUnivar(
            raster=f"single_{test_date}", reference=values, precision=3.1
        )
        values = dict(null_cells=0, min=0, max=18)
        self.assertRasterFitsUnivar(
            raster=f"single_{end_year}_12_31", reference=values, precision=1
        )
        values = dict(mean=0.17)
        self.assertRasterFitsUnivar(
            raster=f"single_{end_year}_12_31", reference=values, precision=0.05
        )

    def test_weather_probabilistic_zero_stddev_lenient(self):
        """Check probabilistic weather with zero stddev and deterministic values"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            weather_coefficient_file=self.weather_file,
            weather_coefficient_stddev_file=self.weather_zero_stddev_file,
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
            overwrite=True,
        )
        test_date = "2021_12_31"
        end_year = end[:4]

        # Values are taken from the deterministic test.

        # Final outputs
        values = dict(null_cells=0, min=0, max=18)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=1)
        values = dict(mean=0.170)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.05)
        values = dict(null_cells=0, min=0, max=100, mean=4.6)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=1.4
        )
        values = dict(null_cells=0, min=0, max=5.46)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.46)
        values = dict(mean=0.25)
        # For single run stddev is zero, but even for many runs the tolerance needs to be
        # almost 100% of the value.
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.25)

        # Time-series outputs
        values = dict(null_cells=0, min=0, max=17.9, mean=0.10)
        self.assertRasterFitsUnivar(
            raster=f"average_{test_date}", reference=values, precision=2.9
        )
        values = dict(null_cells=0, min=0, max=100, mean=3.03)
        self.assertRasterFitsUnivar(
            raster=f"probability_{test_date}", reference=values, precision=0.5
        )
        values = dict(null_cells=0, min=0, max=4.98)
        self.assertRasterFitsUnivar(
            raster=f"stddev_{test_date}", reference=values, precision=0.75
        )
        values = dict(mean=0.14)
        self.assertRasterFitsUnivar(
            raster=f"stddev_{test_date}", reference=values, precision=0.26
        )

        # Single run outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.10)
        self.assertRasterFitsUnivar(
            raster=f"single_{test_date}", reference=values, precision=3.1
        )
        values = dict(null_cells=0, min=0, max=18)
        self.assertRasterFitsUnivar(
            raster=f"single_{end_year}_12_31", reference=values, precision=1
        )
        values = dict(mean=0.17)
        self.assertRasterFitsUnivar(
            raster=f"single_{end_year}_12_31", reference=values, precision=0.05
        )

    def test_weather_probabilistic_stddev_strict(self):
        """Check probabilistic weather"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            weather_coefficient_file=self.weather_file,
            weather_coefficient_stddev_file=self.weather_stddev_file,
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
        )
        test_date = "2021_12_31"
        end_year = end[:4]

        # Values are taken from the deterministic test.

        # Final outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.170)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=4.287)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=5.643, mean=0.127)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

        # Time-series outputs
        values = dict(null_cells=0, min=0, max=18.0, mean=0.101)
        self.assertRasterFitsUnivar(
            raster=f"average_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=100, mean=2.889)
        self.assertRasterFitsUnivar(
            raster=f"probability_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=5.426, mean=0.071)
        self.assertRasterFitsUnivar(
            raster=f"stddev_{test_date}", reference=values, precision=0.001
        )

        # Single run outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.095)
        self.assertRasterFitsUnivar(
            raster=f"single_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=18, mean=0.172)
        self.assertRasterFitsUnivar(
            raster=f"single_{end_year}_12_31", reference=values, precision=0.001
        )

    def test_soil_with_deterministic_weather_strict(self):
        """Check soils with deterministic weather.

        Uses close-to-exact global statistics.
        """
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            weather_coefficient_file=self.weather_file,
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            dispersers_to_soils=0.9,
            soil_survival_steps=100,
            random_seed=1,
            runs=5,
            nprocs=5,
        )
        test_date = "2021_12_31"
        end_year = end[:4]

        # Final outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.040)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=1.120)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=6.000, mean=0.013)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

        # Time-series outputs
        values = dict(null_cells=0, min=0, max=18.0, mean=0.035)
        self.assertRasterFitsUnivar(
            raster=f"average_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=100, mean=1.068)
        self.assertRasterFitsUnivar(
            raster=f"probability_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=5.238, mean=0.011)
        self.assertRasterFitsUnivar(
            raster=f"stddev_{test_date}", reference=values, precision=0.001
        )

        # Single run outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.034)
        self.assertRasterFitsUnivar(
            raster=f"single_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=18, mean=0.040)
        self.assertRasterFitsUnivar(
            raster=f"single_{end_year}_12_31", reference=values, precision=0.001
        )

    def test_nulls_in_input(self):
        """Same as test_outputs() but using inputs with null values."""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host_nulls",
            total_plants="max_host",
            infected="infection_nulls",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")
        end = end[:4]
        self.assertRasterExists(f"average_{end}_12_31")
        self.assertRasterExists(f"probability_{end}_12_31")
        self.assertRasterExists(f"single_{end}_12_31")
        self.assertRasterExists(f"stddev_{end}_12_31")

        values = dict(null_cells=0, min=0, max=18, mean=1.777)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=33.664)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=7.548, mean=0.945)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

    def test_outputs_mortality(self):
        """Check dead output of mortality"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
            flags="m",
            mortality_rate=0.5,
            mortality_time_lag=0,
            mortality_series="dead",
            mortality_frequency="yearly",
        )
        end = end[:4]
        self.assertRasterExists(f"dead_{end}_12_31")

        values = dict(null_cells=0, min=0, max=4.2, mean=0.145)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=8.032)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=16, mean=0.456)
        self.assertRasterFitsUnivar(
            raster=f"dead_{end}_12_31", reference=values, precision=0.001
        )

    def test_outputs_mortality_enabled_rate_0(self):
        """Check with mortality rate 0 (values copied from basic outputs test)"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
            flags="m",
            mortality_rate=0,
            mortality_time_lag=0,
            mortality_series="dead",
            mortality_frequency="yearly",
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")
        end = end[:4]
        self.assertRasterExists(f"average_{end}_12_31")
        self.assertRasterExists(f"probability_{end}_12_31")
        self.assertRasterExists(f"single_{end}_12_31")
        self.assertRasterExists(f"stddev_{end}_12_31")

        ref_float = dict(datatype="DCELL")
        ref_int = dict(datatype="CELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)
        self.assertRasterFitsInfo(raster=f"single_{end}_12_31", reference=ref_int)
        self.assertRasterFitsInfo(raster=f"average_{end}_12_31", reference=ref_float)
        self.assertRasterFitsInfo(
            raster=f"probability_{end}_12_31", reference=ref_float
        )
        self.assertRasterFitsInfo(raster=f"stddev_{end}_12_31", reference=ref_float)

        values = dict(null_cells=0, min=0, max=18, mean=1.777)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=33.664)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=7.547, mean=0.945)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

    def test_outputs_mortality_many_runs(self):
        """Check mortality with many stochastic runs"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=50,
            nprocs=5,
            flags="m",
            mortality_rate=0.5,
            mortality_time_lag=0,
            mortality_series="dead",
            mortality_frequency="yearly",
            overwrite=True,
        )
        end_for_name = end[:4]
        self.assertRasterExists(f"dead_{end_for_name}_12_31")

        # The original reference values were obtained from a run with 100 stochastic
        # runs with seed 1 and 3 non-zero digits were kept. The precision was chosen
        # so that 100 additional runs with 1 stochastic run and different seeds
        # than the original set would each still pass the test.
        # After mortality rounding, the reference values were simply replaced by the
        # new results, but the precision was kept the same. 100 additional runs still
        # pass the test.
        values = dict(null_cells=0, mean=0.143)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.14)
        values = dict(null_cells=0, mean=7.9)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=4.0
        )
        values = dict(null_cells=0, mean=0.456)
        self.assertRasterFitsUnivar(
            raster=f"dead_{end_for_name}_12_31",
            reference=values,
            precision=0.12,
        )

    def test_outputs_mortality_treatment(self):
        """Check mortality together with treatment"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
            mortality_frequency="yearly",
            flags="m",
            mortality_rate=0.5,
            mortality_time_lag=0,
            mortality_series="dead",
            treatments="treatment",
            treatment_date="2020-12-01",
            treatment_length=0,
            treatment_application="ratio_to_all",
        )
        self.assertRasterExists("dead_2019_12_31")
        self.assertRasterExists("dead_2020_12_31")
        self.assertRasterExists("dead_2021_12_31")
        self.assertRasterExists("dead_2022_12_31")

        values = dict(null_cells=0, mean=0.028, sum=1196)
        self.assertRasterFitsUnivar(
            raster="dead_2019_12_31",
            reference=values,
            precision=0.12,
        )
        values = dict(null_cells=0, mean=0.089, sum=2937)
        self.assertRasterFitsUnivar(
            raster="dead_2020_12_31",
            reference=values,
            precision=0.12,
        )
        values = dict(null_cells=0, mean=0.234, sum=5345)
        self.assertRasterFitsUnivar(
            raster="dead_2021_12_31",
            reference=values,
            precision=0.12,
        )
        values = dict(null_cells=0, mean=0.309, sum=8636)
        self.assertRasterFitsUnivar(
            raster="dead_2022_12_31",
            reference=values,
            precision=0.12,
        )

        precision = 0.1
        self.assertRasterFitsUnivar(
            raster="single_2019_12_31", reference=dict(sum=1055), precision=precision
        )
        self.assertRasterFitsUnivar(
            raster="single_2020_12_31", reference=dict(sum=1236), precision=precision
        )
        self.assertRasterFitsUnivar(
            raster="single_2021_12_31", reference=dict(sum=1605), precision=precision
        )
        self.assertRasterFitsUnivar(
            raster="single_2022_12_31", reference=dict(sum=2097), precision=precision
        )

        values = dict(null_cells=0, min=0, max=3.4, mean=0.103)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=5.848)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )

    def test_outputs_mortality_pesticide_treatment(self):
        """Check mortality together with pesticide treatment (all_infected_in_cell)"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
            mortality_frequency="yearly",
            flags="m",
            mortality_rate=0.5,
            mortality_time_lag=0,
            mortality_series="dead",
            treatments="treatment",
            treatment_date="2020-12-01",
            treatment_length=1,
            treatment_application="all_infected_in_cell",
        )
        self.assertRasterExists("dead_2019_12_31")
        self.assertRasterExists("dead_2020_12_31")
        self.assertRasterExists("dead_2021_12_31")
        self.assertRasterExists("dead_2022_12_31")

        values = dict(null_cells=0, mean=0.028)
        self.assertRasterFitsUnivar(
            raster="dead_2019_12_31",
            reference=values,
            precision=0.12,
        )
        values = dict(null_cells=0, mean=0.089)
        self.assertRasterFitsUnivar(
            raster="dead_2020_12_31",
            reference=values,
            precision=0.12,
        )
        values = dict(null_cells=0, mean=0.234)
        self.assertRasterFitsUnivar(
            raster="dead_2021_12_31",
            reference=values,
            precision=0.12,
        )
        values = dict(null_cells=0, mean=0.327)
        self.assertRasterFitsUnivar(
            raster="dead_2022_12_31",
            reference=values,
            precision=0.12,
        )

        values = dict(null_cells=0, min=0, max=3.4, mean=0.111)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=6.160)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )

    def test_sei_treatments_removal_ratio_to_all(self):
        """Check outputs with SEI and treatment_length == 0"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            model_type="SEI",
            latency_period=10,
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
            treatments="treatment",
            treatment_date="2020-12-01",
            treatment_length=0,
            treatment_application="ratio_to_all",
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")
        end = end[:4]
        self.assertRasterExists(f"average_{end}_12_31")
        self.assertRasterExists(f"probability_{end}_12_31")
        self.assertRasterExists(f"single_{end}_12_31")
        self.assertRasterExists(f"stddev_{end}_12_31")

        # Final outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.405)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=10.163)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=6.087, mean=0.343)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

        test_date = "2021_12_31"
        end_year = end[:4]

        # Time-series outputs
        values = dict(null_cells=0, min=0, max=18.0, mean=0.195)
        self.assertRasterFitsUnivar(
            raster=f"average_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=100, mean=5.357)
        self.assertRasterFitsUnivar(
            raster=f"probability_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=5.389, mean=0.147)
        self.assertRasterFitsUnivar(
            raster=f"stddev_{test_date}", reference=values, precision=0.001
        )

        # Single run outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.198)
        self.assertRasterFitsUnivar(
            raster=f"single_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=18, mean=0.405)
        self.assertRasterFitsUnivar(
            raster=f"single_{end_year}_12_31", reference=values, precision=0.001
        )

    def test_sei_treatments_pesticide_ratio_to_all(self):
        """Check outputs with SEI and treatment_length != 0 (pesticide)"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            model_type="SEI",
            latency_period=10,
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
            treatments="treatment",
            treatment_date="2020-12-01",
            treatment_length=500,
            treatment_application="ratio_to_all",
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")
        end = end[:4]
        self.assertRasterExists(f"average_{end}_12_31")
        self.assertRasterExists(f"probability_{end}_12_31")
        self.assertRasterExists(f"single_{end}_12_31")
        self.assertRasterExists(f"stddev_{end}_12_31")

        # Final outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.480)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=11.657)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=6.369, mean=0.388)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

        test_date = "2021_12_31"
        end_year = end[:4]

        # Time-series outputs
        values = dict(null_cells=0, min=0, max=18.0, mean=0.218)
        self.assertRasterFitsUnivar(
            raster=f"average_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=100, mean=6.135)
        self.assertRasterFitsUnivar(
            raster=f"probability_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=5.276, mean=0.164)
        self.assertRasterFitsUnivar(
            raster=f"stddev_{test_date}", reference=values, precision=0.001
        )

        # Single run outputs
        values = dict(null_cells=0, min=0, max=18, mean=0.221)
        self.assertRasterFitsUnivar(
            raster=f"single_{test_date}", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=18, mean=0.484)
        self.assertRasterFitsUnivar(
            raster=f"single_{end_year}_12_31", reference=values, precision=0.001
        )

    def test_outputs_sei_inf(self):
        """Test no change in infected in the first latency period.

        Test that there is no change in infected before reaching end of
        the first latency period.
        """
        start = "2019-01-01"
        end = "2019-01-04"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            output_frequency="daily",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="day",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            model_type="SEI",
            latency_period=10,
            random_seed=1,
            runs=5,
            nprocs=5,
        )
        single = "single_" + end.replace("-", "_")
        self.assertRastersNoDifference(
            reference="infection", actual=single, precision=0
        )

    def test_outputs_sei0(self):
        """Check that outputs of SEI0 have expected values of SI run.

        This is a copy of the basic test_outputs() function except the
        SEI0 addition to parameters.
        """
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            model_type="SEI",
            latency_period=0,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")
        end = end[:4]
        self.assertRasterExists(f"average_{end}_12_31")
        self.assertRasterExists(f"probability_{end}_12_31")
        self.assertRasterExists(f"single_{end}_12_31")
        self.assertRasterExists(f"stddev_{end}_12_31")

        ref_float = dict(datatype="DCELL")
        ref_int = dict(datatype="CELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)
        self.assertRasterFitsInfo(raster=f"single_{end}_12_31", reference=ref_int)
        self.assertRasterFitsInfo(raster=f"average_{end}_12_31", reference=ref_float)
        self.assertRasterFitsInfo(
            raster=f"probability_{end}_12_31", reference=ref_float
        )
        self.assertRasterFitsInfo(raster=f"stddev_{end}_12_31", reference=ref_float)

        values = dict(null_cells=0, min=0, max=18, mean=1.777)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=33.664)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=7.547, mean=0.945)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

    def test_natural_kernel_only(self):
        """Check with only natural kernel (no anthropogenic kernel)"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            random_seed=1,
            runs=5,
            nprocs=5,
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")
        end = end[:4]
        self.assertRasterExists(f"average_{end}_12_31")
        self.assertRasterExists(f"probability_{end}_12_31")
        self.assertRasterExists(f"single_{end}_12_31")
        self.assertRasterExists(f"stddev_{end}_12_31")

        ref_float = dict(datatype="DCELL")
        ref_int = dict(datatype="CELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)
        self.assertRasterFitsInfo(raster=f"single_{end}_12_31", reference=ref_int)
        self.assertRasterFitsInfo(raster=f"average_{end}_12_31", reference=ref_float)
        self.assertRasterFitsInfo(
            raster=f"probability_{end}_12_31", reference=ref_float
        )
        self.assertRasterFitsInfo(raster=f"stddev_{end}_12_31", reference=ref_float)

        values = dict(null_cells=0, min=0, max=18, mean=0.431)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=7.520)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=6.112, mean=0.120)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

    def test_minimal_parameters(self):
        """Check with only minimal set of parameters (optional or with defaults)"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            natural_distance=50,
            random_seed=1,
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")

        ref_float = dict(datatype="DCELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)

        values = dict(null_cells=0, min=0, max=18, mean=2.222)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=44.247)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=0, mean=0)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

    def test_survival_rate_100_percent(self):
        """Check with 100% survival rate

        The values were originally taken from the minimal parameter test.
        """
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            natural_distance=50,
            survival_rate=self.survival_rate_100_percent,
            survival_month=2,
            survival_day=15,
            random_seed=1,
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")

        ref_float = dict(datatype="DCELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)

        values = dict(null_cells=0, min=0, max=18, mean=2.222)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=44.247)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=0, mean=0)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

    def test_survival_rate_0_percent(self):
        """Check with 0% survival rate"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            natural_distance=50,
            survival_rate=self.survival_rate_0_percent,
            survival_month=2,
            survival_day=15,
            random_seed=1,
            runs=3,
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")

        ref_float = dict(datatype="DCELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)

        values = dict(null_cells=0, min=0, max=0, mean=0)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        # Even with multiple runs, stddev should be still zero.
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

    def test_outputs_dispersers(self):
        """Check dead output of mortality"""
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=0.95,
            random_seed=1,
            runs=5,
            nprocs=5,
            dispersers_output="dispersers",
            established_dispersers_output="established_dispersers",
            outside_spores="outside_dispersers",
        )
        self.assertRasterExists("dispersers")
        self.assertRasterExists("established_dispersers")

        values = dict(null_cells=0, min=0, max=15530, mean=522.275)
        self.assertRasterFitsUnivar(
            raster="dispersers", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=129, mean=8.833)
        self.assertRasterFitsUnivar(
            raster="established_dispersers", reference=values, precision=0.001
        )

        self.assertVectorExists("outside_dispersers")
        values = dict(level=2, num_dblinks=0)
        self.assertVectorFitsExtendedInfo(vector="outside_dispersers", reference=values)
        values = dict(points=129936, primitives=129936)
        self.assertVectorFitsTopoInfo(vector="outside_dispersers", reference=values)

    def test_with_and_without_anthropogenic_dispersal_multiple_seeds(self):
        """Check that multiple seeds keep anthropogenic dispersal separate

        Results of running with anthropogenic_dispersal and percent_natural_dispersal=1
        will be the same as without running anthropogenic_dispersal or with
        anthropogenic_dispersal with different seed when multiple seeds are used for
        isolated generators.
        """
        start = "2019-01-01"
        end = "2022-12-31"
        seeds = {
            "disperser_generation": 1,
            "natural_dispersal": 2,
            "anthropogenic_dispersal": 1,
            "establishment": 1,
            "weather": 2,
            "lethal_temperature": 3,
            "movement": 4,
            "overpopulation": 5,
            "survival_rate": 6,
            "soil": 7,
        }
        random_seeds_parameter = (
            ",".join([f"{key}={value}" for key, value in seeds.items()]),
        )
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=1,
            random_seeds=random_seeds_parameter,
            runs=5,
            nprocs=5,
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")
        end_year = end[:4]
        self.assertRasterExists(f"average_{end_year}_12_31")
        self.assertRasterExists(f"probability_{end_year}_12_31")
        self.assertRasterExists(f"single_{end_year}_12_31")
        self.assertRasterExists(f"stddev_{end_year}_12_31")

        # Now we test with a different seed for anthropogenic dispersal.
        # The percent_natural_dispersal is 1, so no anthropogenic dispersal will
        # happen and the result should be the same.
        seeds["anthropogenic_dispersal"] = 2345
        random_seeds_parameter = (
            ",".join([f"{key}={value}" for key, value in seeds.items()]),
        )
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average_2",
            average_series="average_2",
            single_series="single_2",
            stddev="stddev_2",
            stddev_series="stddev_2",
            probability="probability_2",
            probability_series="probability_2",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=1,
            random_seeds=random_seeds_parameter,
            runs=5,
            nprocs=5,
        )

        for reference, actual in [
            ("average", "average_2"),
            (f"single_{end_year}_12_31", f"single_2_{end_year}_12_31"),
            ("stddev", "stddev_2"),
            ("probability", "probability_2"),
        ]:
            self.assertRastersEqual(
                reference,
                actual,
                precision=0.0,
            )

        # Disabling the anthropogenic dispersal completely should still give
        # the same result. The result would not be the same if using only one seed.
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average_3",
            average_series="average_3",
            single_series="single_3",
            stddev="stddev_3",
            stddev_series="stddev_3",
            probability="probability_3",
            probability_series="probability_3",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            random_seeds=random_seeds_parameter,
            runs=5,
            nprocs=5,
        )

        for reference, actual in [
            ("average", "average_3"),
            (f"single_{end_year}_12_31", f"single_3_{end_year}_12_31"),
            ("stddev", "stddev_3"),
            ("probability", "probability_3"),
        ]:
            self.assertRastersEqual(
                reference,
                actual,
                precision=0.0,
            )

        # We test the specific values but they are not important
        # for main purpose of the test.
        ref_float = dict(datatype="DCELL")
        ref_int = dict(datatype="CELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)
        self.assertRasterFitsInfo(raster=f"single_{end_year}_12_31", reference=ref_int)
        self.assertRasterFitsInfo(
            raster=f"average_{end_year}_12_31", reference=ref_float
        )
        self.assertRasterFitsInfo(
            raster=f"probability_{end_year}_12_31", reference=ref_float
        )
        self.assertRasterFitsInfo(
            raster=f"stddev_{end_year}_12_31", reference=ref_float
        )

        values = dict(null_cells=0, min=0, max=18, mean=0.448)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=7.766)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=5.879, mean=0.116)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

    def test_with_and_without_anthropogenic_dispersal_single_seed(self):
        """Check the assumption that single seed does not separate

        Results of running with anthropogenic_dispersal and percent_natural_dispersal=1
        will be different than without running anthropogenic_dispersal.
        """
        start = "2019-01-01"
        end = "2022-12-31"
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average",
            average_series="average",
            single_series="single",
            stddev="stddev",
            stddev_series="stddev",
            probability="probability",
            probability_series="probability",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            anthropogenic_dispersal_kernel="cauchy",
            anthropogenic_distance=1000,
            anthropogenic_direction_strength=0,
            percent_natural_dispersal=1,
            random_seed=1,
            runs=5,
            nprocs=5,
        )
        self.assertRasterExists("average")
        self.assertRasterExists("stddev")
        self.assertRasterExists("probability")
        end_year = end[:4]
        self.assertRasterExists(f"average_{end_year}_12_31")
        self.assertRasterExists(f"probability_{end_year}_12_31")
        self.assertRasterExists(f"single_{end_year}_12_31")
        self.assertRasterExists(f"stddev_{end_year}_12_31")

        # Disabling the anthropogenic dispersal completely should still give
        # the same result.
        self.assertModule(
            "r.pops.spread",
            host="host",
            total_plants="max_host",
            infected="infection",
            average="average_3",
            average_series="average_3",
            single_series="single_3",
            stddev="stddev_3",
            stddev_series="stddev_3",
            probability="probability_3",
            probability_series="probability_3",
            start_date=start,
            end_date=end,
            seasonality=[1, 12],
            step_unit="week",
            step_num_units=1,
            reproductive_rate=1,
            natural_dispersal_kernel="exponential",
            natural_distance=50,
            natural_direction="W",
            natural_direction_strength=3,
            random_seed=1,
            runs=5,
            nprocs=5,
        )

        for reference, actual in [
            ("average", "average_3"),
            (f"single_{end_year}_12_31", f"single_3_{end_year}_12_31"),
        ]:
            # assertRastersNotEqual would be better and we could test all outputs,
            # but it does not exist.
            self.assertRastersDifference(
                reference,
                actual,
                dict(min=-8, max=6),
                precision=8,
            )


if __name__ == "__main__":
    test()
