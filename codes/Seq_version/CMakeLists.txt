#Executable for Misner
set (MISNER_SEQ misner_seq)
set (SOURCES_MISNER_SEQ ${MISNER_SEQ}.cpp)
add_executable(${MISNER_SEQ} ${SOURCES_MISNER_SEQ} ${KADATH_HEADERS} ${CODES_HEADER_FILES})
target_link_libraries(${MISNER_SEQ} kadath ${KADATH_DEPENDENCIES})
target_include_directories(${MISNER_SEQ} PUBLIC ${PROJECT_SOURCE_DIR}/codes/include)

#Executable for Schwarz
set (SCHWARZ_SEQ schwarz_seq)
set (SOURCES_SCHWARZ_SEQ ${SCHWARZ_SEQ}.cpp)
add_executable(${SCHWARZ_SEQ} ${SOURCES_SCHWARZ_SEQ} ${KADATH_HEADERS} ${CODES_HEADER_FILES})
target_link_libraries(${SCHWARZ_SEQ} kadath ${KADATH_DEPENDENCIES})
target_include_directories(${SCHWARZ_SEQ} PUBLIC ${PROJECT_SOURCE_DIR}/codes/include)

#Executable for Kerr
set (KERR_SEQ kerr_seq)
set (SOURCES_KERR_SEQ ${KERR_SEQ}.cpp)
add_executable(${KERR_SEQ} ${SOURCES_KERR_SEQ} ${KADATH_HEADERS} ${CODES_HEADER_FILES})
target_link_libraries(${KERR_SEQ} kadath ${KADATH_DEPENDENCIES})
target_include_directories(${KERR_SEQ} PUBLIC ${PROJECT_SOURCE_DIR}/codes/include)


set_target_properties(${MISNER_SEQ} ${SCHWARZ_SEQ} ${KERR_SEQ} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/sequential)
