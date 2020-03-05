project(Kadath)
cmake_minimum_required(VERSION 2.8)

message("Root directory : ${PROJECT_SOURCE_DIR}")

#set(CMAKE_VERBOSE_MAKEFILE ON)

#Parallel or sequential
option(PAR_VERSION "Parallel version" ON)
option(LINUX "Build for linux" ON)
option(ENABLE_PROFILING "Enable auto-profiling functions." OFF)
option(ENABLE_MOVE_SEMANTIC "Enable move semantic for the Array class" ON)
option(TENSOR_MOVE_SEMANTIC "Enables/disables move ctor and operators fot the Tensor class hierarchy." ON)

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
else(PAR_VERSION)
    #file(COPY ${PROJECT_SOURCE_DIR/src_seq/do_newton.cpp DESTINATION ${PROJECT_SOURCE_DIR/src/System_of_eqs/)
endif(PAR_VERSION)

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


message("Libraries   (from CMakeLists.txt) :")
message(" - threads   :  ${CMAKE_THREAD_LIBS_INIT}")
message(" - GSL       :  ${GSL_LIBRARIES}")
message(" - PGPLOT    :  ${PGPLOT_LIBRARIES}")
message(" - FFTW      :  ${FFTW_LIBRARIES}")
message(" - LAPACK    :  ${LAPACK_LIBRARIES}")
message(" - SCALAPACK :  ${SCALAPACK_LIBRARIES}")

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
if(NOT(PAR_VERSION))
    enable_testing()
endif()
add_subdirectory(tests)