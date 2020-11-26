FABM-NflexPD is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation (https://www.gnu.org/licenses/gpl.html), either version 3 of the License, or (at your option) any later version.
It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
A copy of the license is provided in COPYING.

Copyright 2018 Onur Kerimoglu (kerimoglu.o@gmail.com)

# A few words on NflexPD, the Nutrient- flexible Phytoplankton - Detritus model

This is the model described by Kerimoglu et al., in GMDD, 2020, specific equation numbers in which are referred to in the model code. A previous, non-FABM version was originally published by Smith et al. 2016, Journal of Plankton Research 38, 977–992.

# Obtaining the code and building 
Obtain the FABM-NflexPD code:
    git clone https://github.com/OnurKerimoglu/fabm-nflexpd.git <NFLEXPDDIR>
(Replace `<NFLEXPDDIR>` with the directory with the NFLEXPD code, e.g., ~/nflexpd-git.)

Now download FABM. FABM is a submodule in GOTM, so first download GOTM: 
    git clone https://github.com/gotm-model/code.git <GOTMDIR>
(Replace `<GOTMDIR>` with the directory with the GOTM code, e.g., ~/gotm-git.)
    cd <GOTMDIR>
    git submodule update --init --recursive
(This should download the FABM code under  <GOTMDIR>/extern/fabm)
Check the 'COMPATIBILITY_NOTES' to find out the successfully tested FABM/GOTM versions 
It is advisable to checkout among those compatible before proceeding. 

For letting FABM know about NFLEXPD, edit <GOTMDIR>/extern/fabm/src/CMakeLists.txt and append 
'nflexpd       # Nutrient - flexible Phytoplankton - Detritus' at the end of the DEFAULT_INSTITUTES list. 
Once the model is successfully tested, it is advisable to create a FABM-branch and commit this change. E.g.,
    cd <GOTMDIR>/extern/fabm
    git checkout -b 1.0rc1-nflexpd
    git add src/CMakeLists.txt
    git commit -m 'included nflexpd in CMakeLists'

FABM and NFLEXPD use object-oriented Fortran and therefore require a recent Fortran compiler (gfortran, ifort, etc). Moreover, a platform-independent build system based on [cmake](http://www.cmake.org) is used. Check whether you have that installed: execute `cmake --version` on the command line.

## GOTM + FABM + NFLEXPD

To build GOTM with FABM support, create a build directory, call cmake to generate makefiles, and make to compile and install. For instance:

    mkdir -p ~/build/gotm && cd ~/build/gotm
    cmake <GOTMDIR>/src -DFABM_BASE=<FABMDIR> -DFABM_NFLEXPD_BASE=<NFLEXPDDIR>
    make install

In the above, replace `<GOTMDIR>` with the directory with the GOTM source code, e.g., $HOME/opt/gotm-git if you executed `git clone` in you home directory. Also, replace `<FABMDIR>` with the directory with the FABM code, e.g., <GOTMDIR>/extern/fabm and `<NFLEXPDDIR>` with the directory with the NFLEXPD code, e.g., $HOME/opt/nflexpd-git.

Now you should have a GOTM executable with FABM and NFLEXPD support at `~/local/gotm/bin/gotm`.

It is good practice to keep up to date with the latest code from the NFLEXPD, FABM and GOTM repositories by regularly executing `git pull` in a directory of each repository, and merging these with local branches if necessary.

If either the NFLEXPD, FABM or GOTM source codes change (e.g., because changes you made to the code yourself, or after `git pull`), you will need to recompile. This does NOT require rerunning cmake. Instead, you need to return to the build directory and rerun `make install`. For instance `cd ~/build/gotm && make install`.

## FABM-NFLEXPD 0d

The 0d driver allows you to run FABM models in a "well-mixed box", under arbitrary (time-varying) environmental forcing.

To build the 0d driver, you need to create a directory to build the code in, call `cmake` to generate makefiles, and call `make` to compile and install the FABM library. Usually, the following suffices for this:

    mkdir -p ~/build/fabm-0d && cd ~/build/fabm-0d
    cmake <FABMDIR>/src/drivers/0d -DGOTM_BASE=<GOTMDIR> -DFABM_NFLEXPD_BASE=<NFLEXPDDIR>
    make install

In the above, replace `<NFLEXPDDIR>` with the path to directory with the NFLEXPD source code (e.g., $HOME/opt/nflexpd-git), `<FABMDIR>` with the path to directory with the FABM source code (e.g., ~/fabm-git), and `<GOTMDIR>` with the path to directory with the GOTM source code (e.g., ~/gotm-git). The latter is needed because the 0d driver uses GOTM routines for input, output, time integration, etc.

This will give you an executable at `~/local/fabm/0d/bin/fabm0d`.
