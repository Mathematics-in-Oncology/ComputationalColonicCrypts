## ComputationalColonicCrypts

ComputationalColonicCrypts is a modification and extension of the existing
Chaste implementation on colonic crypt modeling for modeling Lynch syndrome
carcinogenesis.


## Overview
ComputationalColonicCrypts aims to be an implementation of the
computational model presented in Haupt et al (2020), in preparation
for the evolution of colonic crypts in Lynch syndrome carcinogenesis.
The implementation is based on the Chaste software, release 2019_1.
It can be included in the existing Chaste implementation by replacing
some of the existing files and by adding some additional files, all
available here in this repository.

## Documentation

For the Chaste documentation, please read https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides.

### Preliminaries
For running the numerical simulations with Chaste, the following .cpp and .hpp existing
files of Chaste have to be replaced:
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


The following .cpp and .hpp files (if not stated otherwise) have to be added
to the existing Chaste code:
- all mutations except wild-type in
  cell_based/src/cell/properties/mutations
- cell_based/src/population/killers/CryptFissionCellKiller
- cell_based/src/population/killers/HomeostaticCellKiller
- crypt/test/cell/TestColonicCryptSimulation.hpp

### Execute the code
The main file is *crypt/test/cell/TestColonicCryptSimulation.hpp*

To configure and compile the code after replacing and adding the files listed above, go to the build folder and

    cmake ../
    make -j4

Run the code on one processor with

    ctest -V -R TestColonicCryptSimulation

or on multiple processor to simulate multiple crypts in parallel with

    mpirun -n 4 ./run

This will run the code on 4 processors.

### Set parameter values

The following parameter settings can be altered in
*crypt/test/cell/TestColonicCryptSimulation.hpp*:
- the length of the simulation, measured in multiples of the stem cell cycle
- the number of crypts to be simulated sequentially
- the *MLH1*/*MSH2* mutation status
- initial conditions, such as number of cells and possible stem cell mutations
- stem cell mutation rates, exchange rate, death rate
- the number of stem cells
- the crypt fission rate

In *crypt/src/cell/cycle/SimpleWntCellCycleModel.cpp*, we can change
- the *MLH1*/*MSH2* mutation status (this has to be in concordance with the main file)
- effect of mutations on cell differentiation
- mutation rate for TA cells


## System Requirements

### Hardware Requirements
For hardware requirements, we refer to https://chaste.cs.ox.ac.uk/trac/wiki/GettingStarted#SupportedOperatingSystems.

### Software Requirements
For running the simulations, you need to install Chaste, release 2019_1
available on https://www.cs.ox.ac.uk/chaste/download.html.

### OS Requirements
For OS requirements, we refer to https://chaste.cs.ox.ac.uk/chaste/tutorials/release_2019.1/.
We note that the custom code base runs on Ubuntu 18.04 LTS. 

## Installation Guide
The Installation Guide for Chaste is available online https://chaste.cs.ox.ac.uk/chaste/tutorials/release_2019.1/InstallGuides/InstallGuide.html.

Further, the modified and the additional files have to be included in
the Chaste folders as listed above.

## License

Chaste source is released under 3-clause BSD license.

Copyright (c) 2010, The Chancellor, Masters and Scholars of the
University of Oxford. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of the University of Oxford nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
