set(CMAKE_MODULE_PATH $ENV{HOME_KADATH}/Cmake)
if(EXISTS $ENV{HOME_KADATH}/Cmake/CMakeLocal.cmake)
	include ($ENV{HOME_KADATH}/Cmake/CMakeLocal.cmake)	
endif()

option(PAR_VERSION "Parallel version" ON)

file(GLOB_RECURSE HEADERS ${CMAKE_SOURCE_DIR}/include/*.hpp)
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

include_directories($ENV{HOME_KADATH}/include)
include_directories($ENV{HOME_KADATH}/src/Array)
include_directories($ENV{HOME_KADATH}/include/Kadath_point_h)

#Get the correct Kadath library (default == Release)
if (CMAKE_BUILD_TYPE MATCHES Debug)
	set(LIB_KADATH $ENV{HOME_KADATH}/lib/libkadath-debug.a)
else()
	set(LIB_KADATH $ENV{HOME_KADATH}/lib/libkadath.a)
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

#Need scalapack if PAR VERSION
if (PAR_VERSION)
	if (NOT DEFINED SCALAPACK_LIBRARIES)
		message ("Need to give scalapack location in CMakeLocal.cmake")
	endif()
endif()

if(NOT PAR_VERSION)
    if (NOT DEFINED PGPLOT_LIBRARIES)
	    find_package(PGPLOT REQUIRED) #pgplot probably only used in sequential mode
    endif()
endif(NOT PAR_VERSION)

#need to use C++11
add_definitions(-std=c++11)

