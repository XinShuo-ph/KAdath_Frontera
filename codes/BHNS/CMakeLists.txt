cmake_minimum_required(VERSION 2.8)
project(BHNS)

include ($ENV{HOME_KADATH}/Cmake/CMakeExec.cmake)


#Executable 1
set (PROG1 setup)
set (SOURCES1 src/${PROG1}.cpp)
add_executable(${PROG1} ${SOURCES1} ${HEADERS})
target_link_libraries(${PROG1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})

#Executable 2
set (PROG2 solve)
set (SOURCES2 src/${PROG2}.cpp)
add_executable(${PROG2} ${SOURCES2} ${HEADERS})
target_link_libraries(${PROG2} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})
#
#Executable 3
set (PROG3 increase_resolution)
set (SOURCES3 src/${PROG3}.cpp)
add_executable(${PROG3} ${SOURCES3} ${HEADERS})
target_link_libraries(${PROG3} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})
#
#Executable 3
set (PROG4 reader)
set (SOURCES3 src/${PROG4}.cpp)
add_executable(${PROG4} ${SOURCES3} ${HEADERS})
target_link_libraries(${PROG4} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})
#
## overwrite Kadath linking to get shared lib working
#set(LIB_KADATH -Wl,--whole-archive $ENV{HOME_KADATH}/lib/libkadath.a -Wl,--no-whole-archive)
#
#include_directories(/usr/include/python2.7)
#include_directories(/usr/include)

#Library 1
#set (LIB1 kadath_readers)
#set (SOURCES1 src/${LIB1}.cpp)
#add_library(${LIB1} SHARED ${SOURCES1} ${HEADERS})
#target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} -lboost_python)
#set_target_properties(${LIB1} PROPERTIES PREFIX "")
#set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "../lib")
