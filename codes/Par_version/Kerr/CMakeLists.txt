#Executable kerr
set (KERR kerr)
set (SOURCES_KERR ${KERR}.cpp)
add_executable(${KERR} ${SOURCES_KERR} ${HEADERS})
target_link_libraries(${KERR} kadath ${KADATH_DEPENDENCIES})

#Executable reader (for kerr)
if(PGPLOT_FOUND)
    set (READER kerr_reader)
    set (SOURCES_READER reader.cpp)
    add_executable(${READER} ${SOURCES_READER} ${HEADERS})
    target_link_libraries(${READER} kadath ${KADATH_DEPENDENCIES})
    set(KERR_TARGETS ${KERR} ${READER})
else()
    set(KERR_TARGETS ${KERR})
endif()


set_target_properties(${KERR_TARGETS} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/kerr/)