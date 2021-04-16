## ComputationalColonicCrypts

ComputationalColonicCrypts is a modification and extension of the existing
Chaste implementation on colonic crypt modeling for modeling Lynch syndrome
carcinogenesis.


## Overview
ComputationalColonicCrypts aims to be an implementation of the
computational model presented in Haupt et al (2021), in preparation
for the evolution of colonic crypts in Lynch syndrome carcinogenesis.
[A preprint is available on bioRxiv](https://doi.org/10.1101/2020.12.29.424555).
The implementation is based on the Chaste software, release 2019.1.
It can be included in the existing Chaste implementation by replacing
and adding some files, all available here in this repository.

### Citation
If you are using this code in your project, please consider [citing our paper](https://doi.org/10.1101/2020.12.29.424555).

    @article {Haupt2020.12.29.424555,
        author       = {Haupt, Saskia and Gleim, Nils and Ahadova, Aysel and Bl{\"a}ker, Hendrik and von Knebel Doeberitz, Magnus and Kloor, Matthias and Heuveline, Vincent},
        title        = {Computational model investigates the evolution of colonic crypts during Lynch syndrome carcinogenesis},
        elocation-id = {2020.12.29.424555},
        year         = {2020},
        doi          = {10.1101/2020.12.29.424555},
        publisher    = {Cold Spring Harbor Laboratory},
        URL          = {https://www.biorxiv.org/content/early/2020/12/29/2020.12.29.424555},
        eprint       = {https://www.biorxiv.org/content/early/2020/12/29/2020.12.29.424555.full.pdf},
        journal      = {bioRxiv}
    }


## Installation

### Original Chaste software, release 2019.1
In the following we list instructions for installing the Chaste software (release 2019.1) on Ubuntu 18.04 LTS.
For other operating systems, please consult the [Chaste Guides](https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides).

Download the Chaste release 2019.1 from the [GitHub repository](https://github.com/Chaste/Chaste/releases/tag/release_2019.1)

Download the [Chaste Ubuntu package](https://chaste.cs.ox.ac.uk/chaste/tutorials/release_2019.1/InstallGuides/UbuntuPackage.html) via:

    sudo nano /etc/apt/sources.list.d/chaste.list

and, depending on your Ubuntu version, add one of the lines listed on https://chaste.cs.ox.ac.uk/chaste/tutorials/release_2019.1/InstallGuides/UbuntuPackage.html to the `chaste.list` text file. For Ubuntu 18.04 LTS, add the line

    deb http://www.cs.ox.ac.uk/chaste/ubuntu bionic/

Install the Chaste public license key:

    sudo apt-key adv --recv-keys --keyserver hkp://keyserver.ubuntu.com:80 422C4D99

Get the correct dependencies for Chaste:

    sudo apt-get update
    sudo apt-get install chaste-dependencies

Once you have obtained the Chaste source code and installed all dependencies, [you can make a first run of Chaste](https://chaste.cs.ox.ac.uk/chaste/tutorials/release_2019.1/ChasteGuides/CmakeFirstRun.html). Your Chaste source code is in the directory

    /path/to/Chaste_source_code

Create an empty `chaste_build` directory outside the source directory and change into it:

    mkdir /path/to/chaste_build
    cd /path/to/chaste_build

Configure Chaste

    cmake /path/to/Chaste_source_code

If the command runs without any errors, the last lines of the output will be:

    -- Configuring done
    -- Generating done
    -- Build files have been written to: /path/to/chaste_build

Next, compile the Continuous test pack of Chaste, which is the core set of unit tests
that will determine whether everything is working as expected.

    make -j4 Continuous

This command takes about 90 minutes.
If it runs without any errors, one of the last lines of the output will be

    [100%] Built target Continuous

Now run the tests using

    ctest -j4 -L Continuous

which, if successful, prints the following lines towards the end of the output

    100% tests passed, 0 tests failed out of 368

    Label Time Summary:
    ...
    ...

Now, Chaste is installed, and you can go to the next steps to install the Lynch syndrome colonic crypt simulation example.

### Preparing the installation of the Lynch syndrome colonic crypt simulation example with Chaste

Several existing files have to be replaced and further files have to be added to the original Chaste installation.
For this, both folders `cell_based` and `crypt` with subfolders and files
can be downloaded from this repository and moved into the Chaste source code in `/path/to/Chaste_source_code/`.
Note: Do not copy the `run.cpp` in this step.

Merge the folders and replace the already existing files with the files for the Lynch syndrome simulation.
The easiest way to do this is to select both folders `cell_based` and `crypt`,
copy them, go to the Chaste source code, and paste the folders therein.
You will be asked to **merge folders** and to **replace files** which you accept
for all folders and files.

For completeness, we list all files that should be copied or replaced here.

The following existing .cpp and .hpp files of Chaste have to be replaced:

- crypt/src/cell/cycle/SimpleWntCellCycleModel.hpp
- crypt/src/cell/cycle/SimpleWntCellCycleModel.cpp
- cell_based/src/cell/cycle/AbstractPhaseBasedCellCycleModel.hpp
- cell_based/src/cell/cycle/AbstractPhaseBasedCellCycleModel.cpp
- cell_based/src/population/forces/GeneralisedLinearSpringForce.hpp
- cell_based/src/population/forces/GeneralisedLinearSpringForce.cpp
- cell_based/src/simulation/numerical_methods/ForwardEulerNumericalMethod.hpp
- cell_based/src/simulation/numerical_methods/ForwardEulerNumericalMethod.cpp
- cell_based/src/population/AbstractCellPopulation.hpp
- cell_based/src/population/AbstractCellPopulation.cpp
- cell_based/src/population/AbstractCentreBasedCellPopulation.hpp
- cell_based/src/population/AbstractCentreBasedCellPopulation.cpp


The following .cpp and .hpp files have to be added to the original Chaste code:

- cell_based/src/cell/properties/mutations/ApcLOHCellMutationState.hpp
- cell_based/src/cell/properties/mutations/ApcLOHCellMutationState.cpp
- cell_based/src/cell/properties/mutations/ApcTwoHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/BCLOHApcLOHCellMutationState.hpp
- cell_based/src/cell/properties/mutations/BCLOHApcLOHCellMutationState.cpp
- cell_based/src/cell/properties/mutations/BCLOHApcOneHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/BCLOHApcOneHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/BCLOHCellMutationState.hpp
- cell_based/src/cell/properties/mutations/BCLOHCellMutationState.cpp
- cell_based/src/cell/properties/mutations/BCLOHMMR_ApcLOHCellMutationState.hpp
- cell_based/src/cell/properties/mutations/BCLOHMMR_ApcLOHCellMutationState.cpp
- cell_based/src/cell/properties/mutations/BCLOHMMR_ApcOneHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/BCLOHMMR_ApcOneHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/BCLOHMMRTwoHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/BCLOHMMRTwoHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/BCOneHit_ApcLOHCellMutationState.hpp
- cell_based/src/cell/properties/mutations/BCOneHit_ApcLOHCellMutationState.cpp
- cell_based/src/cell/properties/mutations/BCOneHit_ApcOneHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/BCOneHit_ApcOneHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/BetaCateninOneHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/BetaCateninTwoHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/BetaCateninTwoHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/MMR_ApcLOHCellMutationState.hpp
- cell_based/src/cell/properties/mutations/MMR_ApcLOHCellMutationState.cpp
- cell_based/src/cell/properties/mutations/MMR_ApcOneHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/MMR_ApcOneHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/MMR_ApcTwoHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/MMR_ApcTwoHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/MMR_BCOneHit_ApcLOHCellMutationState.hpp
- cell_based/src/cell/properties/mutations/MMR_BCOneHit_ApcLOHCellMutationState.cpp
- cell_based/src/cell/properties/mutations/MMR_BCOneHit_ApcOneHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/MMR_BCOneHit_ApcOneHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/MMR_BCOneHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/MMR_BCOneHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/MMR_BCTwoHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/MMR_BCTwoHitCellMutationState.cpp
- cell_based/src/cell/properties/mutations/MMRTwoHitCellMutationState.hpp
- cell_based/src/cell/properties/mutations/MMRTwoHitCellMutationState.cpp
- cell_based/src/population/killers/CryptFissionCellKiller.hpp
- cell_based/src/population/killers/CryptFissionCellKiller.cpp
- cell_based/src/population/killers/HomeostaticCellKiller.hpp
- cell_based/src/population/killers/HomeostaticCellKiller.cpp
- crypt/test/cell/TestColonicCryptSimulation.hpp

The following .txt file has to be added which contains all the necessary tests

- crypt/test/ContinuousTestPack.txt

### Execute the code
The main file is *crypt/test/cell/TestColonicCryptSimulation.hpp*

If you want to have a new build folder, create an empty `chaste_build_lynch`
directory outside the source directory and change into it:

    mkdir /path/to/chaste_build_lynch
    cd /path/to/chaste_build_lynch

To configure and compile the code after replacing and adding the files listed above, type

    cmake /path/to/Chaste_source_code

As above, if the command runs without any errors, the last lines of the output will be:

    -- Configuring done
    -- Generating done
    -- Build files have been written to: /path/to/chaste_build

Next, compile the Continuous test pack of Chaste

    make -j4 Continuous

Again, as above, this command takes about 90 minutes.
If it runs without any errors, one of the last lines of the output will be

        [100%] Built target Continuous

Run the code on one processor with

    ctest -V -R TestColonicCryptSimulation

### Visualization

For the visualization, we copy the instructions from the [Chaste guide](https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials/RunningCryptSimulationsWithMutations).
Please note that Java has to be installed for this procedure.

To visualize the results, open a new terminal and type

    cd /path/to/Chaste_source_code
    cd anim

To create the java executable

    javac Visualize2dCentreCells.java

Then, to run the visualization, do

    java Visualize2dCentreCells /tmp/$USER/testoutput/Cryptx.y/results_from_time_t

where `x` is the processor the simulation was running on, `y` is the simulated crypt
and `t` is the starting point in time of a stem cell cycle, i.e.\ with the given parameter setting,
a multiple of 1680.

In the results folder `/tmp/$USER/testoutput/Cryptx.y/results_from_time_t`,
there is also a file `cellmutationstates.dat` which tracks the numbers of each
mutation type in the simulation.
These results are tab-separated columns so may be visualized by using
gnuplot, Matlab or similar.

## Modifications

### Run multiple crypt simulations in parallel

If you want to simulate multiple crypts in parallel on multiple processors,
please make sure that MPI is installed.
Place the `run.cpp` file into the `chaste_build_lynch` folder, change into it, compile and run it by

    g++ -o run run.cpp
    mpirun -n 4 ./run

This will run the code on 4 processors.

### Set parameter values

Currently, the parameter setting corresponds to the one in the manuscript.
The parameter setting can be altered.
For this, you can change the following parameters in
*crypt/test/cell/TestColonicCryptSimulation.hpp*:
- the length of the simulation, measured in multiples of the stem cell cycle
- the number of crypts to be simulated sequentially
- the *MLH1*/*MSH2* mutation status
- initial conditions, such as number of cells and possible stem cell mutations
- stem cell mutation rates, exchange rate, death rate
- the number of stem cells
- the crypt fission rate

In *crypt/src/cell/cycle/SimpleWntCellCycleModel.cpp*, you can change
- the *MLH1*/*MSH2* mutation status (this has to be in concordance with the main file)
- effect of mutations on cell differentiation
- mutation rate for TA cells

### Set initial condition
Currently, the active stem cell is initialized with an MMR variant
to simulate the Lynch syndrome scenario.
It can be altered by, e.g., using an APC-mutated stem cell as initial condition.
This is currently a comment in *crypt/test/cell/TestColonicCryptSimulation.hpp*,
starting in line 266

    // APC-mutated stem cell as initial condition.
    ...

These lines of code can be altered to account for other initial mutation states
of the stem cells in a straightforward way.


## System requirements

### Hardware requirements
For hardware requirements, we refer to https://chaste.cs.ox.ac.uk/trac/wiki/GettingStarted#SupportedOperatingSystems.

### Software requirements
For running the simulations, you need to install [Chaste, release 2019.1](https://github.com/Chaste/Chaste/releases/tag/release_2019.1).

### OS requirements
For OS requirements, we refer to https://chaste.cs.ox.ac.uk/chaste/tutorials/release_2019.1/.
The custom code was developed and tested on Ubuntu 18.04 LTS.

## License

Our custom code is released under 3-clause BSD license.

Copyright (c) 2021, Saskia Haupt and Nils Gleim. All rights reserved.

Chaste source is released under 3-clause BSD license.

Copyright (c) 2010, The Chancellor, Masters and Scholars of the
University of Oxford. All rights reserved.

### License text of the 3-clause BSD license

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of the University of Oxford nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
