#Executable for Misner
set (MISNER misner)
set (SOURCES_MISNER ${MISNER}.cpp)
add_executable(${MISNER} ${SOURCES_MISNER} ${HEADERS})
target_link_libraries(${MISNER} ${KLIB} ${KADATH_DEPENDENCIES})

#Executable for Schwarz
set (SCHWARZ schwarz)
set (SOURCES_SCHWARZ ${SCHWARZ}.cpp)
add_executable(${SCHWARZ} ${SOURCES_SCHWARZ} ${HEADERS})
target_link_libraries(${SCHWARZ} ${KLIB} ${KADATH_DEPENDENCIES})

#Executable for Stokes
set (STOKES stokes)
set (SOURCES_STOKES ${STOKES}.cpp)
add_executable(${STOKES} ${SOURCES_STOKES} ${HEADERS})
target_link_libraries(${STOKES} ${KLIB} ${KADATH_DEPENDENCIES})


set_target_properties(${MISNER} ${SCHWARZ} ${STOKES} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
