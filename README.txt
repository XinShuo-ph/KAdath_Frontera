General purpose instructions :

Compiling the library :

Use cmake to generate a build system (better with version 3.1 or higher, but it should work with the later 2.x versions) :
  - Create a build directory where you want the library to be compiled, for exemple at the Kadath directory :
            cd Kadath
            mkdir build
            cd build
  - Run cmake with the Kadath directory as argument and desired compilation options (for an effective release build,
    don't forget the -DCMAKE_BUILD_TYPE=Release option, or Debug for development):
            cmake <options> ..
    In the cases where cmake could not find some required external libraries, one may pass them manually through the
    Kadath/Cmake/CMakeLocal.cmake file (the fftw and scalapack library must usually be provided in that way). Some
    templates are provided for that file in this directory.
  - Run make to compile the library and some samples :
            make
    Use the -j option with a generous (but reasonable) number of process to speed-up the compilation.
    You may also restrain the compilation process to single targets, for exemple "make kadath" will only generate the
    library binary file, "make <code>" with code=kerr, schwarz, bosons, etc. may be used to create the executable of
    the corresponding sample code. "make doc" can be used to generate the Doxygen documentation.

For samples regarding how to compile and link executables using the library, one may use the sources and CMakeLists.txt
of the codes directory as template, or the ones located in the tutorial folder for more independent builds.

The available main cmake build options are the following (the value in parenthesis corresponds to the default settings) :

 -DPAR_VERSION = On/Off (On)            Set to On to build the MPI parallel version of the library, Off for a sequential
                                        version. The MPI parallel version requires ScaLAPack or the Intel Math Kernel
                                        Library (MKL).

 -DMKL_VERSION = On/Off (Off)           Set to On to use the Intel MKL library instead of Scalapack. The cmake detection
                                        script will use two environment variables to find the necessary library
                                        components : INTELROOT and MKLROOT. These variables are set when running the
                                        configuration script provided with the Intel compilers suit and the MKL (
                                        respectively compilervars.sh/csh and mklvars.sh/csh in the bin directory of the
                                        Intel Compilers and MKL installation directory). Use the CMakeLocal.cmake file
                                        in case of dectection issues.

 -DENABLE_GPU_USE = On/Off (Off)        Set to On to build a version of the library that will use a GPU to solve linear
                                        problems when one is available, based on the MAGMA library. If used, Cmake will
                                        need the CUDADIR and MAGMADIR environment variable which provides the paths to
                                        the libraries root directories (a CUDADIR environment variable is needed to
                                        compile MAGMA). Here again, if difficulties are encountered, the CUDA_LIBRARIES
                                        cmake variable may be set directly in the CMakeLocal.cmake file with the
                                        following statement :
                                        set(CUDA_LIBRARIES "-L<path to cuda root dir> -lcublas -lcudart -lcusparse")
                                        Note that in the actual state of the library, GPU are just involved in the
                                        linear solver of the library while, depending on the hardware, the matrix
                                        computation (which is CPU only) is the most expensive part by two order of
                                        magnitude. Thus, one should not use a GPU-enabled Kadath unless he has access
                                        to a node owning both a GPU and a decent amount of CPUs.

 -DUSE_CXX_STANDARD_14 = On/Off (Off)   Allows to downgrade the version of the C++ standard required for the compilation
                                        of the library to 14 instead of the default C++-17. The compilation may
                                        issues some warning if standard 14 is used but the resulting library is stable.

When using the MKL, it is likely that the Intel compiler has to be used for the whole compilation process of the
library. Cmake often doesn't set the MPI compiler consequently, thus one has to provide the path to the MPI wrapper
manually with the -DMPI_CXX_COMPILER and -DMPI_C_COMPILER variables.

Some other, developper-oriented, options are available and may be found in the documentation.


Building an application using Kadath

The easiest way is to use CMake in order to generate a build system for your own applications.
A cmake "header" file with the required Kadath-related variables is provided in the Cmake directory of the library
under the name CMakeExec.cmake. A template CMakeLists.txt file can be found in the tutorial directory, that shows how
to include this CMakeExec.cmake file in your own project. When using this template, one just has to pass the location
of the root directory of the Kadath source files (the Kadath directory if you got the repository from gitlab), as
well as the build directory containing the binaries, using the two respective cmake option :

-DKADATH_SOURCES_DIRECTORY
-DKADATH_BUILD_DIRECTORY

Another option, KADATH_RELATIVE_BUILD_PATH, may be used to specify the build directory using a relative path from the
Kadath root directory.

These can be ommited if you have set the HOME_KADATH environment variable to the Kadath root directory and compiled the
library directly in the same location, using a subdirectory named "build". In the case were the HOME_KADATH environment
variable is set but a custom build directory is used, either KADATH_BUILD_DIRECTORY or KADATH_RELATIVE_BUILD_PATH may
be used respectively to set the absolute path or the relative one.

The CMakeExec.cmake file defines the following variables :
- KADATH_HEADERS : lists the header files of the Kadath library
- KADATH_LIB : the Kadath library binary file
- KADATH_DEPENDENCIES : the library Kadath from which Kadath depends
The following locations are also set as include directory to the project :
- KADATH_SOURCES_DIRECTORY/include
- KADATH_BUILD_DIRECTORY/include
The former is for the library header files, while the later contains a build-specific header with system-specific and
several options that were chosen at the compilation of the library.

Note that the CMakeExec.cmake file will require the same options as the one used for the Kadath library cmake file. One
must then pass the same values that were used during the Kadath build generations for the following options :
-DPAR_VERSION
-DMKL_VERSION
-DENABLE_GPU_USE

At last, use the desired build-type and needed compiler-related options.

Another way to build a custom executable without the CMakeExec.cmake file is to use one of the demonstration
applications

========================================================================================================================
Complementary instructions for MesoPSL

On MesoPSL, the most effective MPI-parallel build is obtained with the Intel compiler and the MKL library. They can
be found in the following locations :
  Intel compilers : /shared/apps/intel/compilers_and_libraries_2019/linux
  MKL library     : /shared/apps/intel/compilers_and_libraries_2019/linux/mkl
The environment configuration script (compilervars.sh and mklvars.sh) may be found in the bin directory of the above
ones, but these two paths may also directly be used as values for the INTELROOT and MKLROOT environment variables.
The MPI compiler wrappers are in the mpi/intel64/bin/ directory of the Intel root dir, the MPI-related cmake options
are then :
 -DMPI_CXX_COMPILER=/shared/apps/intel/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpiicpc
 -DMPI_C_COMPILER=/shared/apps/intel/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpiicc
Here are complete cmake instructions to generate the build system for a release MKL-based parallel-MPI version of Kadath
along with the content of the CMakeLocal.cmake file :

* Pure MPI Case :
    - cmake command line :
            cmake -DCMAKE_BUILD_TYPE=Release -DPAR_VERSION=On -DMKL_VERSION=On -DMPI_CXX_COMPILER=/shared/apps/intel/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpiicpc -DMPI_C_COMPILER=/shared/apps/intel/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpiicc ..
    - CMakeLocal.cmake content :
            set(FFTW_LIBRARIES "/shared/apps/fftw/3.3.8/intelmpi/19.0.0/intel/19.0.0/lib/libfftw3.a")

* MPI / GPU Case :
    - cmake command line :
            cmake -DCMAKE_BUILD_TYPE=Release -DPAR_VERSION=On -DMKL_VERSION=On -DENABLE_GPU_USE=On -DMPI_CXX_COMPILER=/shared/apps/intel/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpiicpc -DMPI_C_COMPILER=/shared/apps/intel/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpiicc ..
    - CMakeLocal.cmake content :
            set(FFTW_LIBRARIES "/shared/apps/fftw/3.3.8/intelmpi/19.0.0/intel/19.0.0/lib/libfftw3.a")
            set(CUDA_LIBRARIES "-L/shared/apps/cuda/10.2/lib64 -lcublas -lcudart -lcusparse")





