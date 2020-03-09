#Executable kerr
set (SCHWARZ_MPI schwarz_mpi)
set (SOURCES_SCHWARZ_MPI ${SCHWARZ_MPI}.cpp)
add_executable(${SCHWARZ_MPI} ${SOURCES_SCHWARZ_MPI} ${HEADERS})
target_link_libraries(${SCHWARZ_MPI} kadath ${KADATH_DEPENDENCIES})