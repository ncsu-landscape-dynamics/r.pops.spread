# r.pops.spread

This is a GRASS GIS module *r.pops.spread* for simulating spread of
pests and pathogens. The module is a GRASS GIS interface to the PoPS
(Pest or Pathogen Spread) model implemented in C++ library maintained
in the [PoPS Core repository](https://github.com/ncsu-landscape-dynamics/pops-core).

The purpose of the *r.pops.spread* module is to provide easy way of
running the model in GRASS GIS environment once you have calibrated
the model for your purposes. You can obtain the calibration from a
colleague or published work or you can calibrate the model manually (in
GRASS GIS) or use the R interface to PoPS called
[rpops](https://github.com/ncsu-landscape-dynamics/rpops) to do that.

The *r.pops.spread* module is available in GRASS GIS Addons repository
and can be installed directly in GRASS GIS through graphical user
interface or using the following command:

```
g.extension r.pops.spread
```

Alternatively, you can obtain latest source code here and install it
from this repository (see below).

Note: Earlier versions of this module were called *r.spread.pest* and
*r.spread.sod*.

## References

Please, in addition to citing or acknowledging this software
(see [CITATION file](CITATION.cff)), cite the following papers:

Ross K. Meentemeyer, Nik J. Cunniffe, Alex R. Cook, Joao A. N. Filipe,
Richard D. Hunter, David M. Rizzo, and Christopher A. Gilligan 2011.
Epidemiological modeling of invasion in heterogeneous landscapes:
spread of sudden oak death in California (1990–2030).
*Ecosphere* 2:art17.
[DOI: 10.1890/ES10-00192.1](https://doi.org/10.1890/ES10-00192.1)

Tonini, Francesco, Douglas Shoemaker, Anna Petrasova, Brendan Harmon,
Vaclav Petras, Richard C. Cobb, Helena Mitasova,
and Ross K. Meentemeyer.
Tangible geospatial modeling for collaborative solutions
to invasive species management.
*Environmental Modelling & Software* 92 (2017): 176-188.
[DOI: 10.1016/j.envsoft.2017.02.020](https://doi.org/10.1016/j.envsoft.2017.02.020)

## Obtaining the latest code

The PoPS library is in a submodule, so use `--recursive` when cloning,
for example:

```
git clone --recursive git@github.com:ncsu-landscape-dynamics/r.pops.spread.git
```

If you have already cloned, you can obtain the submodules using:

```
git submodule update --init
```

## Updating submodule to latest version

To update the submodule, i.e. update submodule's commit used in this
repository, use:

```
git submodule update --remote
```

Note that this change is recorded in the repository. In other words,
the latest commit of a submodule is part of this repository.
The reason for this is that the code in this repository is linked to a
particular commit in the submodule repository (rather than the latest
version). Git works this way to avoid breaking things unexpectedly due
to changes in the submodule repository.

## Updating the code of the submodule

```
cd pops
git checkout master
git add file.hpp
git commit -m "this and that change"
git push
```

```
cd ..
git commit pops -m "update to latest pops commit"
git push
```

## The files

The main.cpp contains the main program to run.

## To run the model

You can use Linux to run the model in the following way.

Open an terminal and install dependencies:

    sudo apt-get install grass-dev

Download the model code as ZIP or using Git:

    git clone ...

Change the current directory to the model directory:

    cd ...

Compile:

    grass --tmp-location XY --exec g.extension module=r.pops.spread url=.

Run:

    grass .../modeling/scenario1 --exec r.pops.spread ...

## Authors

* Francesco Tonini (original R version)
* Zexi Chen (initial C++ version)
* Vaclav Petras (parallelization, GRASS interface, raster handling, ...)
* Anna Petrasova (single species simulation)

See [CHANGELOG.md](CHANGELOG.md) for details about contributions.

## License

This program is free software under the GNU General Public License
(>=v2). Read the file LICENSE for details.
