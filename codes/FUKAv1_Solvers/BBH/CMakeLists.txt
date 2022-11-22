cmake_minimum_required(VERSION 2.8)
project(myproject)

include ($ENV{HOME_KADATH}/Cmake/CMakeExec.cmake)
#Executable 1
set (PROG1 setup)
set (SOURCES1 src/${PROG1}.cpp)
add_executable(${PROG1} ${SOURCES1} ${HEADERS})
target_link_libraries(${PROG1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES})

#Executable 1
set (PROG1 solve)
set (SOURCES1 src/${PROG1}.cpp)
add_executable(${PROG1} ${SOURCES1} ${HEADERS})
target_link_libraries(${PROG1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES})

set (PROG5 increase_resolution)
set (SOURCES5 src/${PROG5}.cpp)
add_executable(${PROG5} ${SOURCES5} ${HEADERS})
target_link_libraries(${PROG5} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES})

set (PROG5 reader)
set (SOURCES5 src/${PROG5}.cpp)
add_executable(${PROG5} ${SOURCES5} ${HEADERS})
target_link_libraries(${PROG5} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES})
