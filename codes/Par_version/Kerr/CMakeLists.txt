#Executable kerr
set (KERR_MPI kerr_mpi)
set (SOURCES_KERR_MPI ${KERR_MPI}.cpp)
add_executable(${KERR_MPI} ${SOURCES_KERR_MPI} ${KADATH_HEADERS} ${CODES_HEADER_FILES})
target_include_directories(${KERR_MPI} PUBLIC ${PROJECT_SOURCE_DIR}/codes/include)
target_link_libraries(${KERR_MPI} kadath ${KADATH_DEPENDENCIES})

#Executable reader (for kerr)
set (KERR_MPI_READER kerr_reader)
if(PGPLOT_FOUND)
    set (SOURCES_KERR_MPI_READER reader.cpp)
    add_executable(${KERR_MPI_READER} ${SOURCES_KERR_MPI_READER} ${KADATH_HEADERS} ${CODES_HEADER_FILES})
    target_link_libraries(${KERR_MPI_READER} kadath ${KADATH_DEPENDENCIES})
    set(KERR_MPI_TARGETS ${KERR_MPI} ${KERR_MPI_READER})
else()
    message("-- PGPLOT unavailable : the ${KERR_MPI_READER} target will be skipped.")
    set(KERR_MPI_TARGETS ${KERR_MPI})
endif()


set_target_properties(${KERR_MPI_TARGETS} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/kerr/)