#Executable kerr
set (SCHWARZ_MPI schwarz_mpi)
set (SOURCES_SCHWARZ_MPI ${SCHWARZ_MPI}.cpp)
add_executable(${SCHWARZ_MPI} ${SOURCES_SCHWARZ_MPI} ${HEADERS})
target_link_libraries(${SCHWARZ_MPI} kadath ${KADATH_DEPENDENCIES})

set_target_properties(${SCHWARZ_MPI} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/schwarz_mpi/)
