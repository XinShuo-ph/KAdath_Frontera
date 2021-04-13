cmake_minimum_required(VERSION 2.8)
project(myproject)

set(KADATH_SOURCES_DIRECTORY $ENV{HOME_KADATH} CACHE PATH "Kadath sources directory")
set(KADATH_RELATIVE_BUILD_PATH build CACHE STRING "Relative path name for Kadath build directory")
set(KADATH_BUILD_DIRECTORY ${KADATH_SOURCES_DIRECTORY}/${KADATH_RELATIVE_BUILD_PATH} CACHE PATH "Kadath directory")

message("Kadath sources directory : ${KADATH_SOURCES_DIRECTORY}")
message("Kadath build directory   : ${KADATH_BUILD_DIRECTORY}")

include(${KADATH_SOURCES_DIRECTORY}/Cmake/CMakeExec.cmake)

#Executable 1
set (PROG1 myprog)
set(PROG1_SOURCES ${PROJECT_SOURCE_DIR}/src/${PROG1}.cpp)
add_executable(${PROG1} ${PROG1_SOURCES} ${KADATH_HEADERS})
target_link_libraries(${PROG1} ${KADATH_LIB} ${KADATH_DEPENDENCIES})
set_target_properties(${PROG1} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

