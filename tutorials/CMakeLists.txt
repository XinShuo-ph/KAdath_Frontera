cmake_minimum_required(VERSION 2.8)
project(myproject)

set(KADATH_SOURCES_DIRECTORY $ENV{HOME_KADATH} CACHE PATH "Kadath sources directory")
set(KADATH_RELATIVE_BUILD_PATH build CACHE STRING "Relative path name for Kadath build directory")
set(KADATH_BUILD_DIRECTORY ${KADATH_SOURCES_DIRECTORY}/${KADATH_BUILD_PATH_NAME} CACHE PATH "Kadath directory")

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

message("Kadath library : ${KADATH_LIB}")

set(STARTER_SOURCES ${PROJECT_SOURCE_DIR}/${STARTER}.cpp)
add_executable(${STARTER} ${STARTER_SOURCES} ${KADATH_HEADERS})
target_link_libraries(${STARTER} ${KADATH_LIB} ${KADATH_DEPENDENCIES})
set_target_properties(${STARTER} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(SCHWARZ_SOURCES ${PROJECT_SOURCE_DIR}/${SCHWARZ}.cpp)
add_executable(${SCHWARZ} ${SCHWARZ_SOURCES} ${KADATH_HEADERS})
target_link_libraries(${SCHWARZ} ${KADATH_LIB} ${KADATH_DEPENDENCIES})
set_target_properties(${SCHWARZ} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

