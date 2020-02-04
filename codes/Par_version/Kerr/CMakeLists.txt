#Executable kerr
set (KERR kerr)
set (SOURCES_KERR ${KERR}.cpp)
add_executable(${KERR} ${SOURCES_KERR} ${HEADERS})
target_link_libraries(${KERR} ${KLIB} ${KADATH_DEPENDENCIES})

#Executable reader (for kerr)
set (READER kerr_reader)
set (SOURCES_READER reader.cpp)
add_executable(${READER} ${SOURCES_READER} ${HEADERS})
target_link_libraries(${READER} ${KLIB} ${KADATH_DEPENDENCIES})


set_target_properties(${KERR} ${READER} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/kerr/)