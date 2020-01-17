#Executable kerr
set (KERR kerr)
set (SOURCES_KERR ${KERR}.cpp)
add_executable(${KERR} ${SOURCES_KERR} ${HEADERS})
target_link_libraries(${KERR} ${KLIB} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})

#Executable reader (for kerr)
set (READER reader)
set (SOURCES_READER ${READER}.cpp)
add_executable(${READER} ${SOURCES_READER} ${HEADERS})
target_link_libraries(${READER} ${KLIB} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})


set_target_properties(${KERR} ${READER} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)