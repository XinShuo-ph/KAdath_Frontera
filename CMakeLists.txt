cmake_minimum_required(VERSION 2.8)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.0)
    cmake_policy(SET CMP0048 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.1)
    cmake_policy(SET CMP0054 NEW)
endif()
project(Kadath VERSION 2.3.0.0)

set(CMAKE_CXX_STANDARD 17)

message("Root directory : ${PROJECT_SOURCE_DIR}")

#The following variable is just for the config.h file to trigger a debug macro flag.
string(TOUPPER ${CMAKE_BUILD_TYPE} STR_DEBUG_MODE)
if(STR_DEBUG_MODE MATCHES DEBUG)
    set(DEBUG_MODE ON)
else()
    set(DEBUG_MODE OFF)
endif()
option(ENABLE_ASSERTS "Enables the assert preprocessor macro (useless in release build mode)" ON)
option(REMOVE_ALL_CHECKS "Remove most of the runtime range or type/tag checks for increase performances." OFF)
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

find_package(Boost)
if(Boost_FOUND)
    set(HAVE_BOOST TRUE)
else()
    set(HAVE_BOOST FALSE)
endif()

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
set(KADATH_HEADERS ${PROJECT_SOURCE_DIR}/include/*.hpp )
message ("Source files" : ${KADATH_SOURCES})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/Kadath_point_h)
#for templates in Array
include_directories(${PROJECT_BINARY_DIR}/include)

if (DEFINED FFTW_INCLUDE_DIR)
    include_directories (${FFTW_INCLUDE_DIR})
endif()


add_library(kadath ${KADATH_SOURCES} include/array_iterator.hpp)
set_target_properties(kadath PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib COMPILE_FLAGS "-fPIC" $
        <$<CONFIG:Debug>:-O0>)

target_include_directories(kadath PUBLIC ${PROJECT_SOURCE_DIR}/include/Kadath_point_h ${PROJECT_BINARY_DIR}/include/)

if (EXISTS ${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
    include(${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
endif()
include (${PROJECT_SOURCE_DIR}/Cmake/CMakeExec.cmake)

configure_file("${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/include/config.h" )

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

find_package(Doxygen)
if(DOXYGEN_FOUND)
    set(DOXYGEN_GENERATE_LATEX YES)
    set(INPUT_FILES_FOR_DOXIGEN ${KADATH_SOURCES} ${KADATH_HEADERS} ${PROJECT_SOURCE_DIR}/doc/kadath.cpp)
    set(DOXYGEN_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/doc)
    set(DOXYGEN_EXCLUDE ${PROJECT_SOURCE_DIR}/include/Kadath_point_h ${PROJECT_SOURCE_DIR}/tests ${PROJECT_SOURCE_DIR}/tutorials ${PROJECT_SOURCE_DIR}/codes)
#    set(DOXYGEN_EXCLUDE_PATTERNS )
    doxygen_add_docs(doc ${INPUT_FILES_FOR_DOXYGEN})
endif()