Cmake command line options for MesoPSL :
cmake -DCMAKE_BUILD_TYPE=Release -DPAR_VERSION=ON -DMKL_VERSION=ON
 -DMPI_CXX_COMPILER=/shared/apps/intel/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpiicpc
 -DMPI_C_COMPILER=/shared/apps/intel/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpiicc
 -DARRAY_MOVE_SEMANTIC=ON -DLINUX=ON -DUSE_MKL_FFTW3_INTERFACE=OFF ..

The last three one are optional (the above passed values are the default one). Note that the FFTW3 MKL interface is
not functional at the moment.

For a sequential build with gcc and gnu libraries :
cmake -DCMAKE_BUILD_TYPE=Release -DPAR_VERSION=OFF

You may run instal_*.sh script for quick configuration and install (for faster compilation, one can edit them
and add the -j option followed by the desired number of compilation threads).