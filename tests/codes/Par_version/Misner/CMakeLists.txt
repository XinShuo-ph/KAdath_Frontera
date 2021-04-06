#Executable Schwarz
set (MISNER_MPI misner_mpi)
set (SOURCES_MISNER_MPI ${MISNER_MPI}.cpp)
add_executable(${MISNER_MPI} ${SOURCES_MISNER_MPI} ${KADATH_HEADERS} ${CODES_HEADER_FILES})
target_link_libraries(${MISNER_MPI} kadath ${KADATH_DEPENDENCIES})
target_include_directories(${MISNER_MPI} PUBLIC ${CMAKE_SOURCE_DIR}/tests/codes/include)

set_target_properties(${MISNER_MPI} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/misner_mpi/)
