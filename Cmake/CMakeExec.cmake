set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/Cmake)
if(EXISTS ${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
	include (${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
endif()

option(PAR_VERSION "Parallel version" ON)
option(MKL_VERSION "MKL Parallel version" OFF)
option(USE_MKL_FFTW3_INTERFACE "Use MKL interface for FFTW3 function (unsupported at the moment)" OFF)

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

if(ENABLE_GPU_USE)
	message("GPU-based computation enabled.")
endif()

file(GLOB_RECURSE HEADERS ${CMAKE_SOURCE_DIR}/include/*.hpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/Kadath_point_h)

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
endif(ENABLE_GPU_USE)

SET(KADATH_DEPENDENCIES  Threads::Threads  ${GSL_LIBRARIES})

if(PAR_VERSION)
	if(MKL_VERSION)
		if(USE_MKL_FFTW3_INTERFACE)
			set(KADATH_DEPENDENCIES ${KADATH_DEPENDENCIES} ${MKL_LIBRARIES})
		else()
			set(KADATH_DEPENDENCIES ${KADATH_DEPENDENCIES} ${FFTW_LIBRARIES} ${MKL_LIBRARIES})
		endif()
	else()
		set(KADATH_DEPENDENCIES ${KADATH_DEPENDENCIES} ${PGPLOT_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})
	endif()
else()
	set(KADATH_DEPENDENCIES ${KADATH_DEPENDENCIES} ${PGPLOT_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

if(ENABLE_GPU_USE)
	set(KADATH_DEPENDENCIES ${MAGMA_LINKER_FLAGS} ${MAGMA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES}
			${CUDA_cusparse_LIBRARY} ${KADATH_DEPENDENCIES})
endif()
#need to use C++17
add_definitions(-std=c++17)

message("===========================================================================")
message("Kadath dependencies : ${KADATH_DEPENDENCIES}")
message("===========================================================================")
