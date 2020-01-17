#cmake_minimum_required(VERSION 2.8)
#project(Kadath_debug)

if (EXISTS ${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
	include(${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
endif()

#Parallel or sequential
option(PAR_VERSION "Parallel version" ON)

if (PAR_VERSION)
	message ("Parallel version")
else(PAR_VERSION)
	message ("Sequential version")
endif(PAR_VERSION)

set(CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR})

if(PAR_VERSION)
#	file(COPY ${PROJECT_SOURCE_DIR/src_par/do_newton.cpp DESTINATION ${PROJECT_SOURCE_DIR/src/System_of_eqs/)
	find_package(MPI REQUIRED)
	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
	message ("MPI Compiler " ${MPI_CXX_COMPILER})
else(PAR_VERSION)
#	file(COPY ${PROJECT_SOURCE_DIR/src_seq/do_newton.cpp DESTINATION ${PROJECT_SOURCE_DIR/src/System_of_eqs/)
endif(PAR_VERSION)

#Get all the src files
file(GLOB_RECURSE KADATH_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp )
message ("Source files" : ${KADATH_SOURCES})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/Kadath_point_h)
#for templates in Array
include_directories(${PROJECT_SOURCE_DIR}/src/Array)
if (DEFINED FFTW_INCLUDE_DIR)
	include_directories (${FFTW_INCLUDE_DIR})
endif()

#need to use C++11
add_definitions(-std=c++11)

add_library(kadath-debug STATIC ${KADATH_SOURCES})
#target_link_libraries(kadath-debug ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})
set_target_properties(kadath-debug PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib COMPILE_FLAGS "-fPIC")

