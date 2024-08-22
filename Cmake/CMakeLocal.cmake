set (PGPLOT_LIBRARIES "/work2/10061/physixin/frontera/miniconda/lib/libpgplot.so")
# set (GSL_LIBRARIES "/opt/apps/intel19/gsl/2.6/lib/libgsl.so")
# set (GSL_LIBRARIES "/work2/10061/physixin/frontera/miniconda/envs/NR/lib/libgsl.so")
set (GSL_LIBRARIES "/opt/apps/intel19/gsl/2.6/lib/libgsl.so")
set (GSL_INCLUDE_DIR  "/opt/apps/intel19/gsl/2.6/include")
include_directories ("/opt/apps/intel19/gsl/2.6/include")
# set (GSL_LIBRARIES "-I/opt/apps/intel19/gsl/2.6/include -L/opt/apps/intel19/gsl/2.6/lib -lgsl -lgslcblas -lm")
# set (SCALAPACK_LIBRARIES "/work2/10061/physixin/frontera/miniconda/envs/NR/lib/libscalapack.so")
# set (SCALAPACK_LIBRARIES "/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so")

set(SCALAPACK_LIBRARIES "-L$ENV{MKLROOT} -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl")
#set (BLACS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblacs-openmpi.so")
# set (FFTW_LIBRARIES "/work2/10061/physixin/frontera/miniconda/envs/NR/lib/libfftw3.so")
set (FFTW_LIBRARIES "/opt/apps/intel19/impi19_0/fftw3/3.3.8/lib/libfftw3.so")

set (FFTW_INCLUDE_DIR "/opt/apps/intel19/impi19_0/fftw3/3.3.8/include")

# set (BLAS_LIBRARIES "/work2/10061/physixin/frontera/miniconda/envs/NR/lib/libblas.so")
# set (BLAS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblas.so")



set (LAPACK_LIBRARIES "/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64/libmkl_lapack95_lp64.a")

# set(SCALAPACK_LIBRARIES "-L$ENV{MKLROOT} -lmkl_lapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl")

