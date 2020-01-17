if(PAR_VERSION)
    add_executable(starter starter_mpi.cpp)
else()
    add_executable(starter starter.cpp)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_link_libraries(starter LINK_PUBLIC kadath-debug)
else()
    target_link_libraries(starter LINK_PUBLIC kadath)
endif()

set_target_properties(starter PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)