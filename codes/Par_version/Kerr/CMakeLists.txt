#Executable kerr
set (KERR kerr)
set (SOURCES_KERR ${KERR}.cpp)
add_executable(${KERR} ${SOURCES_KERR} ${HEADERS})
target_link_libraries(${KERR} kadath ${KADATH_DEPENDENCIES})

#Executable reader (for kerr)
if(PGPLOT_FOUND)
    set (KERR_READER kerr_reader)
    set (SOURCES_KERR_READER reader.cpp)
    add_executable(${KERR_READER} ${SOURCES_KERR_READER} ${HEADERS})
    target_link_libraries(${KERR_READER} kadath ${KADATH_DEPENDENCIES})
    set(KERR_TARGETS ${KERR} ${KERR_READER})
else()
    message("PGPLOT unavailable : the ${KERR_READER} target will be skipped.")
    set(KERR_TARGETS ${KERR})
endif()


set_target_properties(${KERR_TARGETS} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/kerr/)