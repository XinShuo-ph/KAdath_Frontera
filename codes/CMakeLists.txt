cmake_minimum_required(VERSION 2.8)
project(myproject)

include ($ENV{HOME_KADATH}/Cmake/CMakeExec.cmake)

#Executable
add_executable(prog ${SOURCES} ${HEADERS})
target_link_libraries(prog ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})

