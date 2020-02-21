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

#Executable for Stokes
set (SCHWARZ_MULTI_THREAD schwarz-multi-thread)
set (SOURCES_SCHWARZ_MULTI_THREAD ${SCHWARZ_MULTI_THREAD}.cpp)
add_executable(${SCHWARZ_MULTI_THREAD} ${SOURCES_SCHWARZ_MULTI_THREAD} ${HEADERS})
target_link_libraries(${SCHWARZ_MULTI_THREAD} kadath ${KADATH_DEPENDENCIES})


set_target_properties(${MISNER} ${SCHWARZ} ${SCHWARZ_MULTI_THREAD} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
