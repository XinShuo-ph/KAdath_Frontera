#Executable init
set (BOSON_INIT init)
set (SOURCES_BOSON_INIT ${BOSON_INIT}.cpp)
add_executable(${BOSON_INIT} ${SOURCES_BOSON_INIT} ${HEADERS})
target_link_libraries(${BOSON_INIT} kadath ${KADATH_DEPENDENCIES})

#Executable sequence
set (BOSON_SEQUENCE sequence)
set (SOURCES_BOSON_SEQUENCE ${BOSON_SEQUENCE}.cpp)
add_executable(${BOSON_SEQUENCE} ${SOURCES_BOSON_SEQUENCE} ${HEADERS})
target_link_libraries(${BOSON_SEQUENCE} kadath ${KADATH_DEPENDENCIES})

#Executable reader (for boson)
set (BOSON_READER boson_reader)
if(PGPLOT_FOUND)
    set (SOURCES_BOSON_READER reader.cpp)
    add_executable(${BOSON_READER} ${SOURCES_BOSON_READER} ${HEADERS})
    target_link_libraries(${BOSON_READER} kadath ${KADATH_DEPENDENCIES})
    set(BOSON_TARGETS ${BOSON_INIT} ${BOSON_SEQUENCE} ${BOSON_READER})
else()
    message("-- PGPLOT unavailable : the ${BOSON_READER} target will be skipped.")
    set(BOSON_TARGETS ${BOSON_INIT} ${BOSON_SEQUENCE})
endif()


set_target_properties(${BOSON_TARGETS} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/boson/)