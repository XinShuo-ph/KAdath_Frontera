set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/Cmake)
if(EXISTS ${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
	include (${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
endif()

option(PAR_VERSION "Parallel version" ON)
option(MKL_VERSION "Parallel version" OFF)

if (PAR_VERSION)
	message ("Parallel version")
else(PAR_VERSION)
	message ("Sequential version")
endif(PAR_VERSION)

if (MKL_VERSION)
	message ("Version with MKL ")
else(MKL_VERSION)
	message ("Version with scalapack")
endif(MKL_VERSION)

file(GLOB_RECURSE HEADERS ${CMAKE_SOURCE_DIR}/include/*.hpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src/Array)
include_directories(${PROJECT_SOURCE_DIR}/include/Kadath_point_h)

#Get the correct Kadath library (default == Release)
if (CMAKE_BUILD_TYPE MATCHES Debug)
	set(LIB_KADATH ${PROJECT_SOURCE_DIR}/lib/libkadath-debug.a)
else()
	set(LIB_KADATH ${PROJECT_SOURCE_DIR}/lib/libkadath.a)
endif()

#If parallel need to use the PMI wrapper
if(PAR_VERSION)
	find_package (MPI REQUIRED)
	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
	message ("MPI CXX " ${MPI_CXX_COMPILER})
endif(PAR_VERSION)

if (NOT DEFINED GSL_LIBRARIES)
     find_package(GSL REQUIRED) #No FindGSL for cmake 2.8
endif()

if (NOT DEFINED FFTW_LIBRARIES)
     include(FindFFTW) #idem for FindFFTW
endif()

if (NOT DEFINED LAPACK_LIBRARIES)
    include(FindLAPACK) #FindLAPACK does exist
endif()


#include(FindSUNDIALS)

#Need scalapack if PAR VERSION
if (PAR_VERSION)
	if (MKL_VERSION)
		include(FindMKL)
	else(MKL_VERSION)
	if (NOT DEFINED SCALAPACK_LIBRARIES)
		message ("Need to give scalapack location in CMakeLocal.cmake")
	endif()
	endif(MKL_VERSION)
endif(PAR_VERSION)

if(NOT PAR_VERSION)
    if (NOT DEFINED PGPLOT_LIBRARIES)
	    find_package(PGPLOT REQUIRED) #pgplot probably only used in sequential mode
    endif()

    if (NOT DEFINED SUNDIALS_LIBRARIES)
	    include(FindSUNDIALS) #SUNDIAL probably only used in sequential mode	
    endif()
endif(NOT PAR_VERSION)

#need to use C++11
add_definitions(-std=c++11)

