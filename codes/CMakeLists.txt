cmake_minimum_required(VERSION 2.8)
project(myproject)

#The location to Kadath headers and library file can be passed either using the HOME_KADATH environment variable or
#using usual cmake options : cmake -DKADATH_SOURCES_DIRECTORY=<path to kadath root dir> -DKADATH_BUILD_DIRECTORY=<
#path to the library build dir>
set(KADATH_SOURCES_DIRECTORY $ENV{HOME_KADATH} CACHE PATH "Kadath sources directory")
set(KADATH_RELATIVE_BUILD_PATH build CACHE STRING "Relative path name for Kadath build directory")
set(KADATH_BUILD_DIRECTORY ${KADATH_SOURCES_DIRECTORY}/${KADATH_RELATIVE_BUILD_PATH} CACHE PATH "Kadath directory")

message("Kadath sources directory : ${KADATH_SOURCES_DIRECTORY}")
message("Kadath build directory   : ${KADATH_BUILD_DIRECTORY}")

include(${KADATH_SOURCES_DIRECTORY}/Cmake/CMakeExec.cmake)

if(NOT PAR_VERSION)
    set(STARTER starter)
    set(SCHWARZ schwarz)
else()
    set(STARTER starter_mpi)
    set(SCHWARZ schwarz_mpi)
endif()

set(PROG1 myprog1)
set(PROG1_SOURCES ${PROJECT_SOURCE_DIR}/${PROG1}.cpp)
add_executable(${PROG1} ${PROG1_SOURCES} ${KADATH_HEADERS})
target_link_libraries(${PROG1} ${KADATH_LIB} ${KADATH_DEPENDENCIES})

set(PROG2 myprog2)
set(PROG2_SOURCES ${PROJECT_SOURCE_DIR}/${PROG2}.cpp)
add_executable(${PROG2} ${PROG2_SOURCES} ${KADATH_HEADERS})
target_link_libraries(${PROG2} ${KADATH_LIB} ${KADATH_DEPENDENCIES})