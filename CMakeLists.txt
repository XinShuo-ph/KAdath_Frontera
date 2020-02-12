project(Kadath)
cmake_minimum_required(VERSION 2.8)

message("Root directory : ${PROJECT_SOURCE_DIR}")

#set(CMAKE_VERBOSE_MAKEFILE ON)

#Parallel or sequential
option(PAR_VERSION "Parallel version" ON)
option(LINUX "Build for linux" ON)
option(ENABLE_PROFILING "Enable auto-profiling functions." OFF)
option(ENABLE_GPU_USE "Enables the use of GPU(s) if available." OFF)
option(ARRAY_MOVE_SEMANTIC "Enable move semantic for the Array class" ON)

configure_file("${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/include/config.h" )

if (EXISTS ${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
    include(${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
endif()

if (PAR_VERSION)
    message ("Parallel version")
else(PAR_VERSION)
    message ("Sequential version")
endif(PAR_VERSION)

set(CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR})

if(PAR_VERSION)
    #file(COPY ${PROJECT_SOURCE_DIR/src_par/do_newton.cpp DESTINATION ${PROJECT_SOURCE_DIR/src/System_of_eqs/)
    find_package(MPI REQUIRED)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    message ("MPI Compiler " ${MPI_CXX_COMPILER})
endif(PAR_VERSION)

if(ENABLE_GPU_USE)
    find_package(CUDA REQUIRED)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(MAGMA magma)
endif(ENABLE_GPU_USE)

#Get all the src files
file(GLOB_RECURSE KADATH_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp )
message ("Source files" : ${KADATH_SOURCES})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/Kadath_point_h)
#for templates in Array
include_directories(${PROJECT_SOURCE_DIR}/src/Array)
include_directories(${PROJECT_BINARY_DIR}/include)

if (DEFINED FFTW_INCLUDE_DIR)
    include_directories (${FFTW_INCLUDE_DIR})
endif()

#need to use C++17
add_definitions(-std=c++17)


add_library(kadath ${KADATH_SOURCES})
set_target_properties(kadath PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib COMPILE_FLAGS "-fPIC" $
        <$<CONFIG:Debug>:-O0>)

target_include_directories(kadath PUBLIC ${PROJECT_SOURCE_DIR}/include/Kadath_point_h ${PROJECT_BINARY_DIR}/include/)

if (EXISTS ${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
    include(${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
endif()
include (${PROJECT_SOURCE_DIR}/Cmake/CMakeExec.cmake)

if(DEFINED CUDA_INCLUDE_DIR)
    include_directories (${CUDA_INCLUDE_DIRS})
endif()
if(DEFINED MAGMA_INCLUDE_DIRS)
    include_directories(${MAGMA_INCLUDE_DIRS})
endif()

message("Libraries   (from CMakeLists.txt) :")
message(" - GSL       :  ${GSL_LIBRARIES}")
message(" - PGPLOT    :  ${PGPLOT_LIBRARIES}")
message(" - FFTW      :  ${FFTW_LIBRARIES}")
message(" - LAPACK    :  ${LAPACK_LIBRARIES}")
message(" - SCALAPACK :  ${SCALAPACK_LIBRARIES}")
if(ENABLE_GPU_USE)
    message(" - CUDA TOOLKIT : ")
    message("            # CUDA_INCLUDE_DIRS = ${CUDA_INCLUDE_DIRS}")
    message("            # CUDA_LIBRARIES = ${CUDA_LIBRARIES}")
    message("            # CUDA_CUBLAS_LIBRARIES = ${CUDA_CUBLAS_LIBRARIES}")
    message("            # CUDA_cusparse_LIBRARY = ${CUDA_cusparse_LIBRARY}")
    message(" - MAGMA : }")
    message("            # MAGMA_INCLUDE_DIRS = ${MAGMA_INCLUDE_DIRS}")
    message("            # MAGMA_LIBRARIES = ${MAGMA_LIBRARIES}")
endif()


add_subdirectory(tutorials)

if(PAR_VERSION)
    add_subdirectory(codes/Par_version/)
else()
    add_subdirectory(codes/Seq_version)
endif()

#add the tests directory
if(NOT(PAR_VERSION))
    enable_testing()
    add_subdirectory(tests)
endif()