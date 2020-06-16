#Executable kerr
set (SCALA_KERR kerr_sc)
set (SOURCES_SCALA_KERR ${SCALA_KERR}.cpp)
add_executable(${SCALA_KERR} ${SOURCES_SCALA_KERR} ${KADATH_HEADERS} ${CODES_HEADER_FILES})
target_include_directories(${SCALA_KERR} PUBLIC ${PROJECT_SOURCE_DIR}/codes/include)
target_link_libraries(${SCALA_KERR} kadath ${KADATH_DEPENDENCIES})

#Executable Schwarz
#set (SCALA_SCHWARZ schwarz_sc)
#set (SOURCES_SCALA_SCHWARZ ${SCALA_SCHWARZ}.cpp)
#add_executable(${SCALA_SCHWARZ} ${SOURCES_SCALA_SCHWARZ} ${KADATH_HEADERS} ${CODES_HEADER_FILES})
#target_include_directories(${SCALA_SCHWARZ} PUBLIC ${PROJECT_SOURCE_DIR}/codes/include)
#target_link_libraries(${SCALA_SCHWARZ} kadath ${KADATH_DEPENDENCIES})
#
#set (SCALA_BOSON boson_sc)
#set (SOURCES_SCALA_BOSON ${SCALA_BOSON}.cpp)
#add_executable(${SCALA_BOSON} ${SOURCES_SCALA_BOSON} ${KADATH_HEADERS} ${CODES_HEADER_FILES})
#target_include_directories(${SCALA_BOSON} PUBLIC ${PROJECT_SOURCE_DIR}/codes/include)
#target_link_libraries(${SCALA_BOSON} kadath ${KADATH_DEPENDENCIES})

set_target_properties(${SCALA_KERR} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/scalability/)