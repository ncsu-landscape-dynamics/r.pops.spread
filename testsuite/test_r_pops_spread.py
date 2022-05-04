#!/usr/bin/env python3

"""Test of r.pops.spread.

.. moduleauthor:: Anna Petrasova
.. moduleauthor:: Vaclav Petras
"""

import tempfile
from pathlib import Path

from grass.gunittest.case import TestCase
from grass.gunittest.main import test


def items_to_file(items, filename):
    """Save list of items to a file, one item per line"""
    with open(filename, mode="w", encoding="utf-8") as file:
        file.write("\n".join(items))


class TestSpread(TestCase):
    """Tests of r.pops.spread"""
    @classmethod
    def setUpClass(cls):
        """Create input data from the full NC SPM dataset"""
        cls.use_temp_region()
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

        # We can remove the directory only at the end, so we can't use with here
        # for the resources allocated below.
        # pylint: disable=consider-using-with

        cls.tmp_dir = tempfile.TemporaryDirectory()

        years_2019_2022 = 4
        months_in_year = 12
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
            ],
        )

    def tearDown(self):
        """Remove maps after each test method"""
        self.runModule(
            "g.remove",
            flags="f",
            type="raster",
            pattern="average*,single*,stddev*,probability*,dead*",
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
        self.assertRasterExists("average" + "_{}_12_31".format(end))
        self.assertRasterExists("probability" + "_{}_12_31".format(end))
        self.assertRasterExists("single" + "_{}_12_31".format(end))
        self.assertRasterExists("stddev" + "_{}_12_31".format(end))

        ref_float = dict(datatype="DCELL")
        ref_int = dict(datatype="CELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)
        self.assertRasterFitsInfo(
            raster="single" + "_{}_12_31".format(end), reference=ref_int
        )
        self.assertRasterFitsInfo(
            raster="average" + "_{}_12_31".format(end), reference=ref_float
        )
        self.assertRasterFitsInfo(
            raster="probability" + "_{}_12_31".format(end), reference=ref_float
        )
        self.assertRasterFitsInfo(
            raster="stddev" + "_{}_12_31".format(end), reference=ref_float
        )

        values = dict(null_cells=0, min=0, max=18, mean=1.777)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=33.664)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=7.547, mean=0.945)
        self.assertRasterFitsUnivar(raster="stddev", reference=values, precision=0.001)

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
        self.assertRasterExists("average" + "_{}_12_31".format(end))
        self.assertRasterExists("probability" + "_{}_12_31".format(end))
        self.assertRasterExists("single" + "_{}_12_31".format(end))
        self.assertRasterExists("stddev" + "_{}_12_31".format(end))

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
        self.assertRasterExists("dead" + "_{}_12_31".format(end))

        values = dict(null_cells=0, min=0, max=6, mean=0.606)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=24.961)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
        )
        values = dict(null_cells=0, min=0, max=15, mean=0.703)
        self.assertRasterFitsUnivar(
            raster="dead" + "_{}_12_31".format(end), reference=values, precision=0.001
        )

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
        self.assertRasterExists("dead" + "_{}_12_31".format(end_for_name))

        # The reference values were obtained from a run with 100 stochastic runs
        # with seed 1 and 3 non-zero digits were kept. The precision was chosen
        # so that 100 additional runs with 1 stochastic run and different seeds
        # than the original set would each still pass the test.
        values = dict(null_cells=0, mean=0.653)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.14)
        values = dict(null_cells=0, mean=25.9)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=4.0
        )
        values = dict(null_cells=0, mean=0.681)
        self.assertRasterFitsUnivar(
            raster="dead" + "_{}_12_31".format(end_for_name),
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

        values = dict(null_cells=0, min=0, max=6, mean=0.493)
        self.assertRasterFitsUnivar(raster="average", reference=values, precision=0.001)
        values = dict(null_cells=0, min=0, max=100, mean=21.002)
        self.assertRasterFitsUnivar(
            raster="probability", reference=values, precision=0.001
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
        self.assertRasterExists("average" + "_{}_12_31".format(end))
        self.assertRasterExists("probability" + "_{}_12_31".format(end))
        self.assertRasterExists("single" + "_{}_12_31".format(end))
        self.assertRasterExists("stddev" + "_{}_12_31".format(end))

        ref_float = dict(datatype="DCELL")
        ref_int = dict(datatype="CELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)
        self.assertRasterFitsInfo(
            raster="single" + "_{}_12_31".format(end), reference=ref_int
        )
        self.assertRasterFitsInfo(
            raster="average" + "_{}_12_31".format(end), reference=ref_float
        )
        self.assertRasterFitsInfo(
            raster="probability" + "_{}_12_31".format(end), reference=ref_float
        )
        self.assertRasterFitsInfo(
            raster="stddev" + "_{}_12_31".format(end), reference=ref_float
        )

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
        self.assertRasterExists("average" + "_{}_12_31".format(end))
        self.assertRasterExists("probability" + "_{}_12_31".format(end))
        self.assertRasterExists("single" + "_{}_12_31".format(end))
        self.assertRasterExists("stddev" + "_{}_12_31".format(end))

        ref_float = dict(datatype="DCELL")
        ref_int = dict(datatype="CELL")
        self.assertRasterFitsInfo(raster="average", reference=ref_float)
        self.assertRasterFitsInfo(raster="stddev", reference=ref_float)
        self.assertRasterFitsInfo(raster="probability", reference=ref_float)
        self.assertRasterFitsInfo(
            raster="single" + "_{}_12_31".format(end), reference=ref_int
        )
        self.assertRasterFitsInfo(
            raster="average" + "_{}_12_31".format(end), reference=ref_float
        )
        self.assertRasterFitsInfo(
            raster="probability" + "_{}_12_31".format(end), reference=ref_float
        )
        self.assertRasterFitsInfo(
            raster="stddev" + "_{}_12_31".format(end), reference=ref_float
        )

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


if __name__ == "__main__":
    test()
