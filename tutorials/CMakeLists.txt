if(PAR_VERSION)
    add_executable(starter starter_mpi.cpp)
else()
    add_executable(starter starter.cpp)
endif()

target_link_libraries(starter LINK_PUBLIC kadath)

set_target_properties(starter PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)