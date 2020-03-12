#Executable for Misner
set (MISNER misner)
set (SOURCES_MISNER ${MISNER}.cpp)
add_executable(${MISNER} ${SOURCES_MISNER} ${HEADERS})
target_link_libraries(${MISNER} kadath ${KADATH_DEPENDENCIES})

#Executable for Schwarz
set (SCHWARZ schwarz)
set (SOURCES_SCHWARZ ${SCHWARZ}.cpp)
add_executable(${SCHWARZ} ${SOURCES_SCHWARZ} ${HEADERS})
target_link_libraries(${SCHWARZ} kadath ${KADATH_DEPENDENCIES})

#Executable for Schwarze multithread
set (SCHWARZ_MULTI_THREAD schwarz-multi-thread)
set (SOURCES_SCHWARZ_MULTI_THREAD ${SCHWARZ_MULTI_THREAD}.cpp)
add_executable(${SCHWARZ_MULTI_THREAD} ${SOURCES_SCHWARZ_MULTI_THREAD} ${HEADERS})
target_link_libraries(${SCHWARZ_MULTI_THREAD} kadath ${KADATH_DEPENDENCIES})

#Executable for Kerr
set (KERR_SEQUENTIAL kerr_seq)
set (SOURCES_KERR_SEQUENTIAL ${KERR_SEQUENTIAL}.cpp)
add_executable(${KERR_SEQUENTIAL} ${SOURCES_KERR_SEQUENTIAL} ${HEADERS})
target_link_libraries(${KERR_SEQUENTIAL} kadath ${KADATH_DEPENDENCIES})


set_target_properties(${MISNER} ${SCHWARZ} ${SCHWARZ_MULTI_THREAD} ${KERR_SEQUENTIAL} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
