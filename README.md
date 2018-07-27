# r.spread.pest and r.spread.sod

recoding the model to create a c++ version of the SOD-model base on https://github.com/f-tonini/SOD-modeling.
This repository contains the c++ version scripts used to develop a stochastic landscape spread model of forest pathogen *P. ramorum*.

The reference paper: Ross K. Meentemeyer, Nik J. Cunniffe, Alex R. Cook, Joao A. N. Filipe, Richard D. Hunter, David M. Rizzo, and Christopher A. Gilligan 2011. Epidemiological modeling of invasion in heterogeneous landscapes: spread of sudden oak death in California (1990–2030). *Ecosphere* 2:art17. [http://dx.doi.org/10.1890/ES10-00192.1] (http://www.esajournals.org/doi/abs/10.1890/ES10-00192.1) 

## Obtaining the latest code

The PoPSS library is in a submodule, so use `--recursive` when cloning,
for example:

```
git clone --recursive git@github.com:ncsu-landscape-dynamics/r.spread.pest.git
```

If you have already cloned, you can obtain the submodules using:

```
git submodule update --init
```

To update the submodule code together with this repository code use:

```
git pull --recurse-submodules
```

To just update the submodule, use:

```
git submodule update --remote
```

## The files
The main.cpp contains the main program to run.

## To run the model

You can use Linux to run the model in the following way.

Open an terminal and install dependencies:

    sudo apt-get install libgdal-dev libnetcdf-dev

Download the model code as ZIP or using Git:

    git clone ...

Change the current directory to the model directory:

    cd ...

Compile:

    make

Run:

    ./a.out > result.txt

## Authors

* Francesco Tonini (original R version)
* Zexi Chen (initial C++ version)
* Vaclav Petras (parallelization, GRASS interface, raster handling)
* Anna Petrasova (single species simulation)

See [CHANGELOG.md](CHANGELOG.md) for details about contributions.

## License

This program is free software under the GNU General Public License
(>=v2). Read the file LICENSE for details.
