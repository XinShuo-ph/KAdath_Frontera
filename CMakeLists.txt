cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0048 NEW)
project(Kadath VERSION 2.3.0.0)

message("Root directory : ${PROJECT_SOURCE_DIR}")

#set(CMAKE_VERBOSE_MAKEFILE ON)

#User options :
option(PAR_VERSION "Parallel version" ON)
option(ENABLE_GTEST "Enable unit/functional tests using google-test." OFF)
option(HAVE_LINUX "Build for linux" ON)
option(ENABLE_PROFILING "Enable auto-profiling functions." OFF)
option(ENABLE_GPU_USE "Enables the use of GPU(s) if available." OFF)
option(TENSOR_MOVE_SEMANTIC "Enables/disables move ctor and operators fot the Tensor class hierarchy." OFF)
set(POSSIBLE_MEMORY_MAP_TYPES vector unordered_map map flat_map)
set(MEMORY_MAP_TYPE vector CACHE STRING "Memory map container type (possible values are : ${POSSIBLE_MEMORY_MAP_TYPE})")
#set(MEMORY_MAP_TYPE_INDEX 3)
if(MEMORY_MAP_TYPE STREQUAL vector)
    set(MEMORY_MAP_TYPE_INDEX 0)
elseif(MEMORY_MAP_TYPE STREQUAL unordered_map)
    set(MEMORY_MAP_TYPE_INDEX 1)
elseif(MEMORY_MAP_TYPE STREQUAL map)
    set(MEMORY_MAP_TYPE_INDEX 2)
elseif(MEMORY_MAP_TYPE STREQUAL flat_map)
    set(MEMORY_MAP_TYPE_INDEX 3)
endif()
message("Memory map container type : ${MEMORY_MAP_TYPE} (code ${MEMORY_MAP_TYPE_INDEX})")


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
    set(MAGMA_INCLUDE_DIR $ENV{MAGMADIR}/include)
    set(MAGMA_LIB_DIR $ENV{MAGMADIR}/lib)
    set(MAGMA_LINKER_FLAGS "-L${MAGMA_LIB_DIR}")
endif(ENABLE_GPU_USE)

#Get all the src files
file(GLOB_RECURSE KADATH_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp )
message ("Source files" : ${KADATH_SOURCES})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/Kadath_point_h)
#for templates in Array
include_directories(${PROJECT_BINARY_DIR}/include)

if (DEFINED FFTW_INCLUDE_DIR)
    include_directories (${FFTW_INCLUDE_DIR})
endif()

#need to use C++17
add_definitions(-std=c++17)


add_library(kadath ${KADATH_SOURCES} include/array_iterator.hpp)
set_target_properties(kadath PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib COMPILE_FLAGS "-fPIC" $
        <$<CONFIG:Debug>:-O0>)

target_include_directories(kadath PUBLIC ${PROJECT_SOURCE_DIR}/include/Kadath_point_h ${PROJECT_BINARY_DIR}/include/)

if (EXISTS ${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
    include(${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
endif()
include (${PROJECT_SOURCE_DIR}/Cmake/CMakeExec.cmake)

if(DEFINED CUDA_INCLUDE_DIRS)
    include_directories (${CUDA_INCLUDE_DIRS})
endif()
if(DEFINED MAGMA_INCLUDE_DIR)
    include_directories(${MAGMA_INCLUDE_DIR})
endif()

message("Libraries   (from CMakeLists.txt) :")
message(" - threads   :  ${CMAKE_THREAD_LIBS_INIT}")
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
    message(" - MAGMA : ")
    message("            # MAGMA_INCLUDE_DIR = ${MAGMA_INCLUDE_DIR}")
    message("            # MAGMA_LIB_DIR = ${MAGMA_LIB_DIR}")
    message("            # MAGMA_LIBRARIES = ${MAGMA_LIBRARIES}")
endif()

#option to compile with google thread-sanitizer
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -fPIE -pie -g")
#option to compile with google address-sanitizer
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -O1 -fno-omit-frame-pointer -g")

add_subdirectory(tutorials)

if(PAR_VERSION)
    add_subdirectory(codes/Par_version/)
else()
    add_subdirectory(codes/Seq_version)
endif()

#add the tests directory
if(NOT(PAR_VERSION) AND ENABLE_GTEST)
    enable_testing()
endif(NOT(PAR_VERSION) AND ENABLE_GTEST)
add_subdirectory(tests)
