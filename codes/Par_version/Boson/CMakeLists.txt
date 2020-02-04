#Executable init
set (BOSON_INIT init)
set (SOURCES_BOSON_INIT ${BOSON_INIT}.cpp)
add_executable(${BOSON_INIT} ${SOURCES_BOSON_INIT} ${HEADERS})
target_link_libraries(${BOSON_INIT} ${KLIB} ${KADATH_DEPENDENCIES})

#Executable sequence
set (BOSON_SEQUENCE sequence)
set (SOURCES_BOSON_SEQUENCE ${BOSON_SEQUENCE}.cpp)
add_executable(${BOSON_SEQUENCE} ${SOURCES_BOSON_SEQUENCE} ${HEADERS})
target_link_libraries(${BOSON_SEQUENCE} ${KLIB} ${KADATH_DEPENDENCIES})

#Executable reader (for boson)
set (BOSON_READER boson_reader)
set (SOURCES_BOSON_READER reader.cpp)
add_executable(${BOSON_READER} ${SOURCES_BOSON_READER} ${HEADERS})
target_link_libraries(${BOSON_READER} ${KLIB} ${KADATH_DEPENDENCIES})


set_target_properties(${BOSON_INIT} ${BOSON_SEQUENCE} ${BOSON_READER} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/boson/)