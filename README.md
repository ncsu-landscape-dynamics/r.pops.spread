# r.pops.spread

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5160178.svg)](https://doi.org/10.5281/zenodo.5160178)

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

Note: Earlier versions of this module were called *r.spread.pest* and
*r.spread.sod*.

## How to cite

If you use this software or code, please cite the following paper:

* Jones, C., Jones, S., Petrasova, A., Petras, V., Gaydos, D.,
  Skrip, M., Takeuchi, Y., Bigsby, K., and Meentemeyer, R., 2021.
  Iteratively forecasting biological invasions with PoPS and a little help from
  our friends.
  *Frontiers in Ecology and the Environment*
  [DOI: 10.1002/fee.2357](https://doi.org/10.1002/fee.2357)

In case you are using the automatic management feature in rpops or the
steering version of r.pops.spread (from the branch steering), please
cite also:

* Petrasova, A., Gaydos, D.A., Petras, V., Jones, C.M., Mitasova, H. and
  Meentemeyer, R.K., 2020.
  Geospatial simulation steering for adaptive management.
  *Environmental Modelling & Software* 133: 104801.
  [DOI: 10.1016/j.envsoft.2020.104801](https://doi.org/10.1016/j.envsoft.2020.104801)

In addition to citing the above paper, we also encourage you to
reference, link, and/or acknowledge specific version of the software
you are using, for example:

* *We have used r.pops.spread GRASS module version 1.0.2 from
  <https://github.com/ncsu-landscape-dynamics/r.pops.spread>
  [DOI: 10.5281/zenodo.5160179](https://doi.org/10.5281/zenodo.5160179)*.

You can find the DOI for the specific version you are using at
[Zenodo](https://doi.org/10.5281/zenodo.5160178).

## Download

### Download and install

The latest release of the *r.pops.spread* module is available in GRASS GIS Addons repository
and can be installed directly in GRASS GIS through graphical user
interface or using the following command:

```
g.extension r.pops.spread
```

Alternatively, you can obtain latest source code here and install it
from this repository (see below).

### Source code download

Just use Git, but note that the
PoPS Core library is in a submodule, so use `--recursive` when cloning,
for example:

```
git clone --recursive git@github.com:ncsu-landscape-dynamics/r.pops.spread.git
```

If you have already cloned, you can obtain the submodules using:

```
git submodule update --init
```

Note that downloading as ZIP won't include the source code for the submodule,
so downloading as ZIP is not useful for this repo.

## Contributing

Please see the [pops-core](https://github.com/ncsu-landscape-dynamics/pops-core#readme) repository
for contributing best practices and release policies.
Other than that, just open pull requests against this repo.

### Updating submodule to latest version

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

### Updating the code of the submodule

```
cd pops-core
git switch -c new-feature-branch
git add file.hpp
git commit -m "this and that change"
git push
```

Then create a PR. After the PR is merged, then update the submodule and commit:

```
cd ..
git submodule update --remote
git commit pops-core -m "update to latest pops commit"
git push
```

## The files

The `main.cpp` file contains the main program to run.
The model is in `pops-core/include/pops-core` directory.

## To run the model

You can use Linux to run the model in the following way.

Open an terminal and install dependencies:

    sudo apt-get install grass grass-dev

Download this repo using Git (see above):

    git clone ...

Change the current directory to the model directory:

    cd ...

Compile:

    grass --tmp-location XY --exec g.extension module=r.pops.spread url=.

Run (assuming you checked how to create a GRASS GIS mapset with our data):

    grass .../modeling/scenario1 --exec r.pops.spread ...

## Authors and contributors

### Authors

_(alphabetical order)_

* Chris Jones
* Margaret Lawrimore
* Vaclav Petras
* Anna Petrasova

### Previous contributors

_(alphabetical order)_

* Zexi Chen
* Devon Gaydos
* Francesco Tonini

See Git commit history, GitHub insights, or CHANGELOG.md file (if present)
for details about contributions.

## License

Permission to use, copy, modify, and distribute this software and its documentation
under the terms of the GNU General Public License version 2 or higher is hereby
granted. No representations are made about the suitability of this software for any
purpose. It is provided "as is" without express or implied warranty.
See the GNU General Public License for more details.
