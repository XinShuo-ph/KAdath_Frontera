if(EXISTS ${KADATH_SOURCES_DIRECTORY}/Cmake/CMakeLocal.cmake)
	include (${KADATH_SOURCES_DIRECTORY}/Cmake/CMakeLocal.cmake)
endif()
set(CMAKE_MODULE_PATH ${KADATH_SOURCES_DIRECTORY}/Cmake)

option(PAR_VERSION "Parallel version" ON)
option(MKL_VERSION "MKL Parallel version" OFF)
option(HAVE_LINUX "Build for linux" ON)
option(ENABLE_GPU_USE "Enables the use of GPU(s) if available." OFF)
option(USE_CXX_STANDARD_14 "Compile the library with C++-14 instead of 17" OFF)
if(USE_CXX_STANDARD_14)
	set(CPP_STANDARD "14")
else()
	set(CPP_STANDARD "17")
endif()

if(NOT KADATH_LIBRARY_CMAKE_BUILD)
	message ("  ")
	message("Kadath-specific CMake options summary : ")
	if(PAR_VERSION)
		message(" - PAR_VERSION ....................: ON  - building MPI parallel version")
	else()
		message(" - PAR_VERSION ....................: OFF - building sequential version")
	endif()
	message(" - HAVE_LINUX .....................: ${HAVE_LINUX}")
	if(ENABLE_GPU_USE)
		message(" - ENABLE_GPU_USE .................: ON  - build hybrid Magma-based GPU / MPI version")
	else()
		message(" - ENABLE_GPU_USE .................: OFF - do not build for GPU usage")
	endif()
	message(" - CPP_STANDARD ...................: ${CPP_STANDARD}")
else(NOT KADATH_LIBRARY_CMAKE_BUILD)
	message("MPI-related Kadath options :")
endif(NOT KADATH_LIBRARY_CMAKE_BUILD)
if(MKL_VERSION)
	message(" - MKL_VERSION ....................: ON  - build will links with intel MKL MPI parallel linear solvers")
else()
	message(" - MKL_VERSION ....................: OFF - link with SCALAPACK to solve linear systems")
endif()
message("   ")

file(GLOB_RECURSE KADATH_HEADERS ${CMAKE_SOURCE_DIR}/include/*.hpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE})

include_directories(${KADATH_SOURCES_DIRECTORY}/include)
include_directories(${KADATH_SOURCES_DIRECTORY}/include/Kadath_point_h)
include_directories(${KADATH_BUILD_DIRECTORY}/include)

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

if(NOT DEFINED CMAKE_THREAD_LIBS_INIT)
	find_package(Threads REQUIRED)
endif()

#include(FindSUNDIALS)

#Need scalapack if PAR VERSION
if (PAR_VERSION)
	if (MKL_VERSION)
		include(FindMKL)
	else(MKL_VERSION)
		if (NOT DEFINED SCALAPACK_LIBRARIES)
			message ("WARNING : MPI parallel version without MKL and SCALAPACK not found.")
			message("           To solve this issue, you may consider passing the SCALAPCK ")
			message("           library location whithin the CMake/CMakeLocal.cmake file.")
		endif()
	endif(MKL_VERSION)
endif(PAR_VERSION)


if (NOT DEFINED PGPLOT_LIBRARIES)
	if(NOT PAR_VERSION)
		find_package(PGPLOT REQUIRED) #pgplot probably only used in sequential mode
	else()
		find_package(PGPLOT) #in this case, the pgplot lib is not directly used by the lib
	endif()
endif()

if(NOT PAR_VERSION)
    if (NOT DEFINED SUNDIALS_LIBRARIES)
	    include(FindSUNDIALS) #SUNDIAL probably only used in sequential mode	
    endif()
endif(NOT PAR_VERSION)

if(ENABLE_GPU_USE)
	if(NOT DEFINED CUDA_LIBRARIES)
		find_package(CUDA REQUIRED)
	endif()
	if(NOT DEFINED MAGMA_LIBRARIES)
		set(MAGMA_INCLUDE_DIR $ENV{MAGMADIR}/include)
		set(MAGMA_LIB_DIR $ENV{MAGMADIR}/lib)
		set(MAGMA_LINKER_FLAGS -L${MAGMA_LIB_DIR})
		set(MAGMA_LIBRARIES "${MAGMA_LINKER_FLAGS} -lmagma")
	endif()
	include_directories(${MAGMA_INCLUDE_DIR})
	include_directories(${CUDA_INCLUDE_DIRS})
endif(ENABLE_GPU_USE)

SET(KADATH_DEPENDENCIES  ${GSL_LIBRARIES})

if(PAR_VERSION)
	if(MKL_VERSION)
		set(KADATH_DEPENDENCIES ${KADATH_DEPENDENCIES} ${FFTW_LIBRARIES} ${MKL_LIBRARIES})
	else()
		set(KADATH_DEPENDENCIES ${KADATH_DEPENDENCIES} ${PGPLOT_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES})
	endif()
	set(KADATH_DEPENDENCIES ${KADATH_DEPENDENCIES} ${SCALAPACK_LIBRARIES})
else()
	set(KADATH_DEPENDENCIES ${KADATH_DEPENDENCIES} ${PGPLOT_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

if(ENABLE_GPU_USE)
	set(KADATH_DEPENDENCIES ${MAGMA_LINKER_FLAGS} ${MAGMA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES}
			${CUDA_cusparse_LIBRARY} ${KADATH_DEPENDENCIES})
endif()
#need to use C++17 (or at least 14)
set(CMAKE_CXX_STANDARD ${CPP_STANDARD})
if(CMAKE_CXX_STANDARD GREATER_EQUAL 17)
	set(USE_CXX_STANDARD_17_OR_HIGHER ON)
else()
	set(USE_CXX_STANDARD_17_OR_HIGHER OFF)
endif()

set(KADATH_LIB ${KADATH_BUILD_DIRECTORY}/lib/libkadath.a)
#for backward compatibility
set(LIB_KADATH ${KADATH_LIB})
set(HEADERS ${KADATH_HEADERS})
#message("  ")
#message("===========================================================================")
#message("Kadath dependencies : ${KADATH_DEPENDENCIES}")
#message("===========================================================================")
