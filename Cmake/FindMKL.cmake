# - Find Intel MKL
# Find the MKL libraries
#
# Options:
#
#   MKL_STATAIC       :   use static linking
#   MKL_MULTI_THREADED:   use multi-threading
#   MKL_SDL           :   Single Dynamic Library interface
#   USE_MKL_FFTW3_INTERFACE : Enable MKL FFT through FFTW3 interface
#
# This module defines the following variables:
#
#   MKL_FOUND            : True if MKL_INCLUDE_DIR are found
#   MKL_INCLUDE_DIR      : where to find mkl.h, etc.
#   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES        : the library to link against.
#   MKL_LIB_WITHOUT_FFT  : the library without MKL FFT.


include(FindPackageHandleStandardArgs)

set(INTEL_ROOT $ENV{INTELROOT} CACHE PATH "Folder contains intel libs")
set(MKL_ROOT $ENV{MKLROOT} CACHE PATH "Folder contains MKL")

# Find include dir
find_path(MKL_INCLUDE_DIR mkl.h
    PATHS ${MKL_ROOT}/include)

# Find include directory
#  There is no include folder under linux
if(WIN32)
    find_path(INTEL_INCLUDE_DIR omp.h
        PATHS ${INTEL_ROOT}/include)
    set(MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR} ${INTEL_INCLUDE_DIR})
endif()

# Find libraries

# Handle suffix
set(_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if(WIN32)
    if(MKL_STATAIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES _dll.lib)
    endif()
else()
    if(MKL_STATAIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
    endif()
endif()


# MKL is composed by four layers: Interface, Threading, Computational and RTL
if(LINUX)
    set(MKL_SEARCH_PATH ${MKL_ROOT}/lib/intel64_lin)
else()
    set(MKL_SEARCH_PATH ${MKL_ROOT}/lib/intel64)
endif()

if(MKL_SDL)
    find_library(MKL_LIBRARY mkl_rt
        PATHS ${MKL_SEARCH_PATH})

    set(MKL_MINIMAL_LIBRARY ${MKL_LIBRARY})
else()
    ######################### Interface layer #######################
    if(WIN32)
        set(MKL_INTERFACE_LIBNAME mkl_intel_c)
    else()
        set(MKL_INTERFACE_LIBNAME mkl_intel_lp64)
    endif()

    find_library(MKL_INTERFACE_LIBRARY ${MKL_INTERFACE_LIBNAME}
        PATHS ${MKL_SEARCH_PATH})

    ######################## Threading layer ########################
    if(MKL_MULTI_THREADED)
        set(MKL_THREADING_LIBNAME mkl_intel_thread)
    else()
        set(MKL_THREADING_LIBNAME mkl_sequential)
    endif()

    find_library(MKL_THREADING_LIBRARY ${MKL_THREADING_LIBNAME}
        PATHS ${MKL_SEARCH_PATH})

    ####################### Computational layer #####################
    find_library(MKL_CORE_LIBRARY mkl_core
        PATHS ${MKL_SEARCH_PATH})
    if(USE_MKL_FFTW3_INTERFACE)
        find_library(MKL_FFT_LIBRARY mkl_cdft_core
            PATHS ${MKL_SEARCH_PATH})
    endif()
    find_library(MKL_SCALAPACK_LIBRARY mkl_scalapack_lp64
        PATHS ${MKL_SEARCH_PATH})
    find_library(MKL_BLAS_LIBRARY mkl_blacs_intelmpi_lp64
        PATHS ${MKL_SEARCH_PATH})

    ############################ RTL layer ##########################
    if(WIN32)
        set(MKL_RTL_LIBNAME libiomp5md)
    else()
        set(MKL_RTL_LIBNAME iomp5)
    endif()
    if(LINUX)
        set(INTEL_RTL_SEARCH_PATH ${INTEL_ROOT}/lib/intel64_lin)
    else()
        set(INTEL_RTL_SEARCH_PATH ${INTEL_ROOT}/lib/intel64)
    endif()
    find_library(MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME}
        PATHS ${INTEL_RTL_SEARCH_PATH})

    if(USE_MKL_FFTW3_INTERFACE)
        set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_FFT_LIBRARY} ${MKL_BLAS_LIBRARY} ${MKL_SCALAPACK_LIBRARY} ${MKL_RTL_LIBRARY})
    else()
        set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_BLAS_LIBRARY} ${MKL_SCALAPACK_LIBRARY} ${MKL_RTL_LIBRARY})
    endif()
    set(MKL_MINIMAL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY})

    message("MKL_Libraries found :")
    message("   MKL_INTERFACE_LIBRARY=${MKL_INTERFACE_LIBRARY}")
    message("   MKL_THREADING_LIBRARY= ${MKL_THREADING_LIBRARY}")
    message("   MKL_CORE_LIBRARY= ${MKL_CORE_LIBRARY}")
    if(USE_MKL_FFTW3_INTERFACE)
        message("   MKL_FFT_LIBRARY= ${MKL_FFT_LIBRARY}")
    endif()
    message("   MKL_SCALAPACK_LIBRARY =  ${MKL_SCALAPACK_LIBRARY}")
    message("   MKL_RTL_LIBRARY= ${MKL_RTL_LIBRARY}")
    message("   MKL_BLAS_LIBRARY= ${MKL_BLAS_LIBRARY}")
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

find_package_handle_standard_args(MKL DEFAULT_MSG
    MKL_INCLUDE_DIR MKL_LIBRARY MKL_MINIMAL_LIBRARY)

if(MKL_FOUND)
    set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
    set(MKL_LIBRARIES ${MKL_LIBRARY})
    set(MKL_MINIMAL_LIBRARIES ${MKL_LIBRARY})
endif()
