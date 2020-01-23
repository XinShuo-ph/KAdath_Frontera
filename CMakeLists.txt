project(Kadath)
cmake_minimum_required(VERSION 2.8)

message("Root directory : ${PROJECT_SOURCE_DIR}")

#set(CMAKE_VERBOSE_MAKEFILE ON)

#Parallel or sequential
option(PAR_VERSION "Parallel version" ON)
option(ARRAY_MOVE_SEMANTIC "Enable move semantic for the Array class" ON)

configure_file("${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/include/config.h" )


if(CMAKE_BUILD_TYPE MATCHES Debug)
    include(build_debug/CMakeLists.txt)
    target_include_directories(kadath-debug PUBLIC ${PROJECT_SOURCE_DIR}/include/Kadath_point_h ${PROJECT_BINARY_DIR}/include/)
    set(KLIB kadath-debug)
else()
    include(build_release/CMakeLists.txt)
    target_include_directories(kadath PUBLIC ${PROJECT_SOURCE_DIR}/include/Kadath_point_h ${PROJECT_BINARY_DIR}/include/)
    set(KLIB kadath)
endif()

if (EXISTS ${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
    include(${PROJECT_SOURCE_DIR}/Cmake/CMakeLocal.cmake)
endif()
include (${PROJECT_SOURCE_DIR}/Cmake/CMakeExec.cmake)


message("Libraries   (from CMakeLists.txt) :")
message(" - GSL       :  ${GSL_LIBRARIES}")
message(" - PGPLOT    :  ${PGPLOT_LIBRARIES}")
message(" - FFTW      :  ${FFTW_LIBRARIES}")
message(" - LAPACK    :  ${LAPACK_LIBRARIES}")
message(" - SCALAPACK :  ${SCALAPACK_LIBRARIES}")


add_subdirectory(tutorials)

if(PAR_VERSION)
    add_subdirectory(codes/Par_version/)
else()
    add_subdirectory(codes/Seq_version)
endif()

#add the tests directory
enable_testing()
if(NOT(PAR_VERSION))
    enable_testing()
    add_subdirectory(tests)
endif()