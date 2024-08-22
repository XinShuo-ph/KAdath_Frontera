\mainpage FUKA Reference
# Frankfurt University/Kadath Initial Data branch
#### Author(s)    : L. Jens Papenfort, Samuel D. Tootle, Philippe Grandclément

## Overview - Updated
Included are the Frankfurt initial data solvers and utilities based on the Kadath
  spectral solver library.  The original solvers written by the aformentioned authors
  (hereafter denoted as FUKAv1) are located in ./codes/FUKAv1/[BH, NS, BHNS, BNS, BBH] respectively.
  The solvers in the next version release, denoted as FUKAv2, can be found in ./codes/FUKAv2/[BH, NS, BHNS, BNS, BBH] respectively.
	Both FUKAv1 (v1) and FUKAv2 (v2) includes support for polytropic equations of state as well as tabulated EOS
  in the standard LORENE format.  Examples and additional details can be found in the [eos](https://bitbucket.org/fukaws/fuka/src/fuka/eos/) directory.

## FUKAv2.2 Release Notes:
With this release, user demanded shells around compact objects have been removed.  In its place is a significantly more robust solution
which determines an optimal configuration of spherical grids around each compact object in binary spaces. 
  
## FUKA Maintainer(s):  

Samuel D. Tootle - tootle@itp.uni-frankfurt.de

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
in `src/Utilities/Exporters`.  These exporters allow one to compile an interface code 
for an evolution framework along with the Kadath static library located in `$HOME_KADATH/lib`, in 
order to export data based on input grid points.  

# 2. Modifications from base Kadath

In addition to the solving routines included within `$HOME_KADATH/codes`, we also note the major overall modifications 
and additions that differ from base Kadath.  
1.  This branch includes memory optimizations that inspired portions of the optimization (now main) branch  
2.  Modification/addition of numerical spaces for the BH, BBH, BNS, and BHNS  
3.  Addition of an equation of state infrastructure utilizing Margherita standalone to handle
tabulated and polytropic EOS - see [include/EOS](https://bitbucket.org/fukaws/fuka/src/fuka/include/EOS)  
4.  Addition of the Configurator framework to enable extensibility of solvers by managing controls,
stages, and key variables - see [include/Configurator](https://bitbucket.org/fukaws/fuka/src/fuka/include/Configurator)  
5.  Addition of exporters for all the previously mentioned ID types - see [src/Utilities/Exporters](https://bitbucket.org/fukaws/fuka/src/fuka/src/Utilities/Exporters)

**Note: as of summer 2021, the FUKA solvers are based on the deprecated branch of Kadath.  Given the optimizations and changes made
within the FUKA branch conflict with those implimented in the `master` branch (previously the `optimized` branch), a considerable
level of effort is required to merge FUKA with the new `master` branch as well as test to see which optimizations provide better results.
Currently, there is no timeline for when this will be done.**

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
2. Create a build directory and enter it
3. Invoke `cmake (options) ..`
  - where `..` denotes the location where the CMakeList.txt file is
  - The important cmake options are the following (the value in parentheses corresponds to the default settings) :
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
one must specify them manually through the [CMakeLocal.cmake](https://bitbucket.org/fukaws/fuka/src/fuka/Cmake/CMakeLocal.cmake) file 
(the `fftw` and `scalapack` libraries must usually be provided in this way). 
Some working [CMakeLocal.cmake](https://bitbucket.org/fukaws/fuka/src/fuka/Cmake/CMakeLocal.cmake) files 
are provided for HPC systems in Germany as well as examples for personal computers.

Once cmake has been successfully invoked, use make -j $KAD_NUMC to start the compilation.

## Compiling the library with the compile script

A script called [compile](https://bitbucket.org/fukaws/fuka/src/fuka/build_release/compile) is also provided that can be used to facilitate the installation process. So long as the above environment variables are set and the libraries are found, no additional input is necessary.
Run the compile script within the `build_release` directory using 

`. compile` 

in order to build the library.  

## Compiling FUKA solvers
The above mentioned [compile script](https://bitbucket.org/fukaws/fuka/src/fuka/build_release/compile) has been added as a symbolic link to the FUKAv1 and FUKAv2 solver directories for convenience to compile the individual solvers.

# 5. Dependencies

1. C++ compiler - gcc is recommended
    - Must support c++17 standards
    - Must support `<filesystem>`
2. Cmake
3. git
4. FFTW3
5. GSL
6. scaLAPACK
7. MPI
8. Boost
9. Boost::python (to compile [PythonTools](https://bitbucket.org/fukaws/fuka/src/fuka/codes/PythonTools/))



# For build on Frontera

set the Kadath env

```bash
export HOME_KADATH=$WORK/Kadath
export KAD_CC=gcc
export KAD_CXX=g++
export KAD_NUMC=32
```

compilers:
```shell
(base) c202-030[clx](607)$ which gcc
/opt/apps/gcc/9.1.0/bin/gcc
(base) c202-030[clx](608)$ which g++
/opt/apps/gcc/9.1.0/bin/g++
(base) c202-030[clx](609)$ which gfortran
/opt/apps/gcc/9.1.0/bin/gfortran
(base) c202-030[clx](616)$ which cmake
/opt/apps/cmake/3.24.2/bin/cmake
(base) c202-030[clx](610)$ which mpirun
/opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpirun
```


modules:

```shell
(base) c209-005[clx](551)$ module list

Currently Loaded Modules:
  1) git/2.24.1      3) cmake/3.24.2    5) xalt/2.10.34   7) gcc/9.1.0   9) python3/3.8.2  11) fftw3/3.3.8
  2) autotools/1.2   4) hwloc/1.11.12   6) TACC           8) gsl/2.6    10) impi/19.0.9    12) mkl/19.0.5
```

download pgplot with conda:

```
conda install conda-forge::pgplot
```


check the paths in `Kadath/Cmake/CMakeLocal.cmake`
```Makefile
set (PGPLOT_LIBRARIES "/work2/10061/physixin/frontera/miniconda/lib/libpgplot.so")
set (GSL_LIBRARIES "/opt/apps/intel19/gsl/2.6/lib/libgsl.so")
set (GSL_INCLUDE_DIR  "/opt/apps/intel19/gsl/2.6/include")
include_directories ("/opt/apps/intel19/gsl/2.6/include")
set(SCALAPACK_LIBRARIES "-L$ENV{MKLROOT} -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl")
set (FFTW_LIBRARIES "/opt/apps/intel19/impi19_0/fftw3/3.3.8/lib/libfftw3.so")
set (FFTW_INCLUDE_DIR "/opt/apps/intel19/impi19_0/fftw3/3.3.8/include")
set (LAPACK_LIBRARIES "/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64/libmkl_lapack95_lp64.a")
```

check the cmake output
```
-- The C compiler identification is GNU 9.1.0
-- The CXX compiler identification is GNU 9.1.0
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working C compiler: /opt/apps/gcc/9.1.0/bin/gcc - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /opt/apps/gcc/9.1.0/bin/g++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
Parallel version
Version with scalapack
-- Found MPI_C: /opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/lib/release/libmpi.so (found version "3.1") 
-- Found MPI_CXX: /opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/lib/libmpicxx.so (found version "3.1") 
-- Found MPI: TRUE (found version "3.1")  
MPI CXX /opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpigxx
CMAKE_CXX_FLAGS_DEBUG is -g -O0 -DNDEBUG
CMAKE_CXX_FLAGS_RELEASE is -g -O3 -std=c++17
-- Configuring done
-- Generating done
```