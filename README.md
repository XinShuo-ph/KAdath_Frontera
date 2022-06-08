# Frankfurt University/Kadath Initial Data branch
#### Author(s)    : L. Jens Papenfort, Samuel D. Tootle, Philippe Grandclément

## Overview - Updated
Included are the Frankfurt initial data solvers and utilities based on the Kadath
  spectral solver library.  The original solvers written by the aformentioned authors
  (hereafter denoted as FUKAv1) are located in ./codes/FUKAv1/[BH, NS, BHNS, BNS, BBH] respectively.
  The solvers in the next version release, denoted as FUKAv2, can be found in ./codes/FUKAv2/[BH, NS, BHNS, BNS, BBH] respectively.
	Both FUKAv1 (v1) and FUKAv2 (v2) includes support for polytropic equations of state as well as tabulated EOS in
  in the standard LORENE format.  Examples and additional details can be found in the [eos](https://bitbucket.org/fukaws/fuka/src/fukav2/eos/) directory.

## FUKAv2 Notes: 
The release of FUKAv2 is considerable step forward in reliable generation of extremal spinning, asymmetric
binary initial data using the KADATH spectral software.  The v2 solvers aim to maximize convergence by using super-imposed
isolated solutions to setup the initial guess for binary ID.  Additionally, v2 aims to automate the generation of ID by minimizing
the workflow for the user.  Finally, v2 allows quite for considerable flexibility in setting up the config file to make ID generation 
as efficient as possible within the KADATH framework including reusing previously solved implicit isolated solutions.  For more details,
please see the documentation in the [FUKAv2](https://bitbucket.org/fukaws/fuka/src/fukav2/codes/FUKAv2/).

  - FUKAv1 Specific: There have been some core changes to various utilities used when constructing ID that have required some refactoring of the v1 solvers.  Spot testing has been done to ensure these codes function as originally intended - however - these codes, overall, remain unchanged.
  
## FUKA Maintainer(s):  

Samuel D. Tootle - tootle@itp.uni-frankfurt.de (primary),  

L. Jens Papenfort - papenfort@itp.uni-frankfurt.de  

## KADATH Maintainer:
Philippe Grandclément - philippe.grandclement@obspm.fr

License      : GPLv3+ for all other code  

## REQUIRED CITATIONS:

1) L. Jens Papenfort, Samuel D. Tootle, Philippe Grandclément, Elias R. Most, Luciano Rezzolla: https://arxiv.org/abs/2103.09911  
  
2) Philippe Grandclément, http://dx.doi.org/10.1016/j.jcp.2010.01.005  

# 1. Purpose

This collection of ID solvers aims at delivering consistent initial data (ID)
solutions to the eXtended Conformal Thin-Sandwich (XCTS) formulation of
Einstein's field equations for a variety of compact object configurations.
  
As each solver has their own specific nuances and considerations, we have included
a README in each solver directory to provide a basis for getting started with the 
respective solver.  
  
Additionally, each initial data has a respective exporter which can be seen 
in src/Utilities/Exporters.  These exporters allow one to compile an interface code 
for an evolution toolkit along with the Kadath static library located in ./lib, in 
order to export data based on input grid points.  

# 2. Modifications from base Kadath

In addition to the solving routines included within ./codes, we also note the major overall modifications 
and additions that differ from base Kadath.  
1.  This branch includes memory optimizations that inspired portions of the optimization branch  
2.  Modification/addition of numerical spaces for the BH, BBH, BNS, and BHNS  
3.  Addition of an equation of state infrastructure utilizing Margherita standalone to handle
tabulated and polytropic EOS - see [include/EOS](https://bitbucket.org/fukaws/fuka/src/fukav2/include/EOS)  
4.  Addition of the Configurator framework to enable extensibility of solvers by managing controls,
stages, and key variables - see [include/Configurator](https://bitbucket.org/fukaws/fuka/src/fukav2/include/Configurator)  
5.  Addition of exporters for all the previously mentioned ID types - see [src/Utilities/Exporters](https://bitbucket.org/fukaws/fuka/src/fukav2src/Utilities/Exporters)

Note: as of summer 2021, the FUKA solvers are based on the deprecated branch of Kadath.  Given the optimizations and changes made
within the FUKA branch conflict with those implimented in the `master` branch (previously the `optimized` branch), a considerable
level of effort is required to merge FUKA with the new `master` branch as well as test to see which optimizations provide better results.
Currently, there is no timeline for when this will be done.

# 3. Public Thorns for use with the Einstein Toolkit

The following workspace includes the FUKA ID respository (including versioned branches) 
as well as available thorns for use with the Einstein Toolkit in order to import FUKA ID:
https://bitbucket.org/fukaws/



# 4. Compiling KADATH

## Get the Sources

`git clone git@bitbucket.org:fukaws/fuka.git`

## Set the following environment variables based on your setup in your ~/.bashrc file

- HOME_KADATH - e.g., export HOME_KADATH=/home/user/lib/fuka/
- KAD_CC - e.g. gcc
- KAD_CXX - e.g. g++
- KAD_NUMC - number of parallel compiling jobs CMake can run, e.g. 7

## Build Process
1. Go to build_release.
1. Create a build directory and enter it
1. Invoke `cmake (options) ..`
  - where `..` denotes the location where the CMakeList.txt file is
  - The available main cmake options are the following (the value in parentheses corresponds to the default settings) :
    - `-DPAR_VERSION = On/Off (On)`
      - Set to On to build the MPI parallel version of the library. The initial data codes within this branch are only designed for use with MPI.
    - `-DCMAKE_BUILD_TYPE = Release/Debug`
      - Specifies the build type which results in different compiler options being used
    - `-DMPI_CXX_COMPILER = <compiler pile>`
      - Path to the MPI C++ wrapper (when not automatically detected by cmake)
    - `-DMPI_C_COMPILER  = <compiler pile>`
      - Path to the MPI C wrapper (when not automatically detected by cmake)

Example using GNU+mpi compilers:  
    `cmake -DCMAKE_BUILD_TYPE=Release -DPAR_VERSION=On -DMPI_CXX_COMPILER=mpic++ -DMPI_C_COMPILER=mpicc ..`

In most HPC systems, `cmake` will likely not find the dependency libraries that the user may intend.  Therefore,
one must specify them manually through the [CMakeLocal.cmake](https://bitbucket.org/fukaws/fuka/src/fukav2/Cmake/CMakeLocal.cmake) file 
(the `fftw` and `scalapack` libraries must usually be provided in this way). 
Some working [CMakeLocal.cmake](https://bitbucket.org/fukaws/fuka/src/fukav2/Cmake/CMakeLocal.cmake) files are provided for HPC systems in Germany as well
as examples for personal computers.

Once cmake has been successfully invoked, use make -j $KAD_NUMC to start the compilation.

## Compiling the library with the compile script

A script called [compile](https://bitbucket.org/fukaws/fuka/src/fukav2/build_release/compile) is also provided that can be used to facilitate the installation process. So long as the above environment variables are set and the libraries are found, no additional input is necessary.
Run the compile script within the `build_release` directory using 

`. compile` 

in order to build the library.  

## Compiling FUKA solvers
The above mentioned [compile script](https://bitbucket.org/fukaws/fuka/src/fukav2/build_release/compile) has been added as a symbolic link to the FUKAv1 and FUKAv2 solver directories for convenience in compiling the individual solvers.

# 5. Dependencies

1. C++ compiler - gcc is recommended
    - Must support c++17 standards
    - Must support `<filesystem>`
1. Cmake
1. git
1. FFTW3
1. GSL
1. scaLAPACK
1. MPI