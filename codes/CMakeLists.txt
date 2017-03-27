cmake_minimum_required(VERSION 2.8)
project(myprog)

include ($ENV{HOME_KADATH}/Cmake/CMakeExec.cmake)

add_executable(myprog ${SOURCES} ${HEADERS})
target_link_libraries(myprog ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})
