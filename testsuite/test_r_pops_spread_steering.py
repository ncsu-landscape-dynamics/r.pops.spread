#!/usr/bin/env python3

"""Test of the computational steering interface of r.pops.spread.

The steering code is reached only when the ip_address and port options are
set, so none of it is covered by test_r_pops_spread.py. These tests drive the
tool over its steering TCP protocol the same way a steering client does.

.. moduleauthor:: Anna Petrasova
"""

import os
import socket
import subprocess
import tempfile
import time

import grass.script as gs
from grass.gunittest.case import TestCase
from grass.gunittest.main import test

# The tool connects to the client, so the test is the server. Waiting for a
# tool which crashed or hung must not block the test run forever.
ACCEPT_TIMEOUT = 60
RECEIVE_TIMEOUT = 300
EXIT_TIMEOUT = 120


class SteeringSession:
    """Drives one run of r.pops.spread over the steering protocol.

    Used as a context manager so that the socket is closed and the tool is
    reaped even when an assertion fails in the middle of a session.
    """

    def __init__(self, **parameters):
        self.parameters = parameters
        self.received = ""
        self._server = None
        self._connection = None
        self._process = None
        self._stderr = None

    def __enter__(self):
        self._server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        # Port 0 lets the OS pick a free port, so parallel tests cannot clash.
        self._server.bind(("127.0.0.1", 0))
        self._server.listen(1)
        self._server.settimeout(ACCEPT_TIMEOUT)
        port = self._server.getsockname()[1]

        arguments = dict(self.parameters, ip_address="127.0.0.1", port=port)
        flags = arguments.pop("flags", "")
        command = ["r.pops.spread", "--overwrite", "--quiet"]
        command += [f"-{flags}"] if flags else []
        command += [
            f"{key}={','.join(str(i) for i in value) if isinstance(value, (list, tuple)) else value}"
            for key, value in arguments.items()
        ]
        # Keep stderr so that an aborted tool can be reported with its message
        # instead of just a closed connection.
        self._stderr = tempfile.TemporaryFile()
        self._process = subprocess.Popen(
            command, stdout=subprocess.DEVNULL, stderr=self._stderr
        )
        self._server.settimeout(1)
        for _ in range(ACCEPT_TIMEOUT):
            try:
                self._connection, _ = self._server.accept()
                break
            except socket.timeout:
                # A tool which rejected the parameters is never going to
                # connect, so report its error instead of waiting.
                if self._process.poll() is not None:
                    raise AssertionError(
                        f"Tool did not connect: {self._tool_error()}"
                    ) from None
        else:
            raise AssertionError("Timed out waiting for the tool to connect")
        self._connection.settimeout(RECEIVE_TIMEOUT)
        return self

    def __exit__(self, *unused):
        if self._connection is not None:
            self._connection.close()
        if self._server is not None:
            self._server.close()
        if self._process is not None and self._process.poll() is None:
            self._process.kill()
            self._process.wait()
        if self._stderr is not None:
            self._stderr.close()

    def _tool_error(self):
        """Message describing how the tool ended, for use in failures"""
        try:
            # The tool may still be on its way out, wait so that the exit code
            # can be reported instead of None.
            code = self._process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            code = None
        self._stderr.seek(0)
        errors = self._stderr.read().decode("utf-8", "replace").strip()
        tail = os.linesep.join(errors.splitlines()[-5:])
        return f"tool exited with {code}{os.linesep}{tail}"

    def send(self, message):
        """Send one or more commands, e.g. 'goto:0;cmd:play;'

        Commands sent together are queued by the tool in order, so no waiting
        between them is needed.
        """
        self._connection.sendall(message.encode("utf-8"))

    def wait_for(self, token, count=1):
        """Block until token has been received count times in total"""
        while self.received.count(token) < count:
            try:
                data = self._connection.recv(4096)
            except socket.timeout:
                raise AssertionError(
                    f"Timed out waiting for {token!r}: {self._tool_error()}"
                ) from None
            if not data:
                raise AssertionError(
                    f"Connection closed waiting for {token!r}: {self._tool_error()}"
                )
            self.received += data.decode("utf-8", "replace")

    def outputs(self):
        """Names of rasters the tool reported through the output messages"""
        names = []
        for part in self.received.split("output:")[1:]:
            names.append(part.split("|")[0])
        return names

    def finish(self):
        """Stop the tool and check that it ended cleanly"""
        self.send("cmd:stop;")
        try:
            code = self._process.wait(timeout=EXIT_TIMEOUT)
        except subprocess.TimeoutExpired:
            raise AssertionError("Tool did not exit after stop command") from None
        if code != 0:
            raise AssertionError(f"Steering run failed: {self._tool_error()}")


class TestSteering(TestCase):
    """Tests of the steering interface of r.pops.spread"""

    # Kept small on purpose, the steering runs are the slow part of the suite.
    common = dict(
        host="host",
        total_plants="max_host",
        infected="infection",
        start_date="2019-01-01",
        end_date="2022-12-31",
        seasonality=[1, 12],
        step_unit="week",
        step_num_units=1,
        output_frequency="yearly",
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
        runs=2,
        nprocs=2,
    )
    mortality = dict(
        flags="m",
        single_series="steer_single",
        mortality_rate=0.5,
        mortality_time_lag=0,
        mortality_series="dead",
        mortality_frequency="yearly",
    )

    @classmethod
    def setUpClass(cls):
        """Create input data from the full NC SPM dataset"""
        cls.use_temp_region()
        cls.temp_dir = (
            tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        )
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
            "v.to.rast", input="railroads", output="infection_", use="val", value=1
        )
        cls.runModule("r.null", map="infection_", null=0)
        cls.runModule("r.mapcalc", expression="infection = if(host > 0, infection_, 0)")
        cls.runModule("r.mapcalc", expression="max_host = 100")
        cls.runModule(
            "r.circle",
            output="raw_infected_patch",
            coordinates=[639300, 220900],
            max=100,
            flags="b",
        )
        cls.runModule(
            "r.mapcalc", expression="infected_patch = min(raw_infected_patch, host)"
        )
        cls.runModule(
            "g.region", n="n-800", s="s+800", e="e-800", w="w+800",
            align="lsat7_2002_30",
        )
        cls.runModule("r.mapcalc", expression="quarantine = 1")
        cls.runModule("g.region", raster="lsat7_2002_30", res=85.5, flags="a")

    @classmethod
    def tearDownClass(cls):
        cls.temp_dir.cleanup()
        cls.del_temp_region()
        cls.runModule(
            "g.remove",
            flags="f",
            type="raster",
            name=[
                "ndvi", "host", "infection_", "infection", "max_host",
                "raw_infected_patch", "infected_patch", "quarantine",
            ],
        )

    def tearDown(self):
        self.runModule(
            "g.remove", flags="f", type="raster", pattern="steer_*,plain_*,dead_*"
        )

    def test_playthrough_matches_plain_run(self):
        """Steering without stepping back gives the same result as no steering"""
        self.assertModule("r.pops.spread", average="plain_average", **self.common)
        with SteeringSession(average="steer_average", **self.common) as session:
            session.send("cmd:play;")
            session.wait_for("info:last:")
            session.finish()
        self.assertRastersNoDifference(
            actual="steer_average", reference="plain_average", precision=0
        )

    def test_goto_start_and_replay_with_mortality(self):
        """Going back to the start and replaying keeps the state consistent

        Mortality cohorts accumulate over the simulation. If they are not
        restored together with the infected hosts, PoPS Core reports more
        dead hosts than there are infected ones and the tool aborts.
        """
        parameters = dict(self.common, **self.mortality)
        with SteeringSession(average="steer_average", **parameters) as session:
            session.send("cmd:play;")
            session.wait_for("info:last:")
            session.send("goto:0;cmd:play;")
            session.wait_for("info:last:", count=2)
            session.finish()

    def test_step_back_and_replay_with_mortality(self):
        """Stepping back one checkpoint and replaying keeps the state consistent"""
        parameters = dict(self.common, **self.mortality)
        with SteeringSession(average="steer_average", **parameters) as session:
            session.send("cmd:play;")
            session.wait_for("info:last:")
            session.send("cmd:stepb;cmd:play;")
            session.wait_for("info:last:", count=2)
            session.finish()

    def test_reported_outputs_exist(self):
        """Rasters announced through the protocol are actually created"""
        with SteeringSession(
            single_series="steer_single", average_series="steer_average", **self.common
        ) as session:
            session.send("cmd:play;")
            session.wait_for("info:last:")
            session.finish()
        names = session.outputs()
        self.assertTrue(names, msg="No output was reported over the protocol")
        existing = gs.list_strings(type="raster", mapset=".")
        for name in names:
            self.assertIn(f"{name}@{gs.gisenv()['MAPSET']}", existing)


    def quarantine_rows(self, path):
        """Data rows currently in the quarantine CSV"""
        if not os.path.exists(path):
            return []
        with open(path, encoding="utf-8") as file:
            return [
                line
                for line in file.read().splitlines()
                if line and not line.startswith("step,")
            ]

    def wait_for_quarantine_rows(self, path, count, timeout=10):
        """Wait until the quarantine file holds the given number of data rows

        The file is rewritten as the simulation advances, so a read can land
        while it is being written.
        """
        deadline = time.time() + timeout
        rows = self.quarantine_rows(path)
        while len(rows) != count and time.time() < deadline:
            time.sleep(0.05)
            rows = self.quarantine_rows(path)
        return rows

    def test_quarantine_is_readable_while_stepping(self):
        """Quarantine results are available for each step as it is computed

        The file is needed while stepping through the simulation, not only
        once the whole simulation ends.
        """
        path = os.path.join(self.temp_dir.name, "quarantine.csv")
        parameters = dict(
            self.common,
            infected="infected_patch",
            single_series="steer_single",
            quarantine="quarantine",
            quarantine_output=path,
        )
        with SteeringSession(**parameters) as session:
            for step in range(1, 5):
                session.send("cmd:stepf;")
                session.wait_for("output:", count=step)
                self.assertEqual(
                    len(self.wait_for_quarantine_rows(path, step)),
                    step,
                    msg=f"Expected {step} quarantine row(s) after {step} step(s)",
                )
            session.finish()

    def test_quarantine_does_not_report_steps_rolled_back(self):
        """Going back in time drops quarantine results of the discarded steps"""
        path = os.path.join(self.temp_dir.name, "quarantine_rollback.csv")
        parameters = dict(
            self.common,
            infected="infected_patch",
            single_series="steer_single",
            quarantine="quarantine",
            quarantine_output=path,
        )
        with SteeringSession(**parameters) as session:
            session.send("cmd:play;")
            session.wait_for("info:last:")
            full = len(self.quarantine_rows(path))
            self.assertGreater(full, 1, msg="Expected several quarantine steps")
            session.send("goto:1;")
            session.finish()
        self.assertEqual(
            len(self.quarantine_rows(path)),
            1,
            msg="Quarantine file still reports steps which were rolled back",
        )


if __name__ == "__main__":
    test()
