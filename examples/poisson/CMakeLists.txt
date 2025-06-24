# Check if we have the necessary source files
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
    # Add the poisson example executable
    add_executable(crossflow_poisson
        main.cpp
    )

    # Set include directories
    target_include_directories(crossflow_poisson 
        PRIVATE 
            ${CMAKE_SOURCE_DIR}/src
    )
    
    # Link against the crossflow library and VTK if available
    target_link_libraries(crossflow_poisson
        PRIVATE
            crossflow
    )
    
    if(VTK_FOUND)
        target_link_libraries(crossflow_poisson PRIVATE ${VTK_LIBS})
        target_compile_definitions(crossflow_poisson PRIVATE -DHAVE_VTK)
    endif()

    # Copy the executable to the binary directory
    set_target_properties(crossflow_poisson PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
    
    # Create output directory during build
    add_custom_command(TARGET crossflow_poisson POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin/out
        COMMENT "Creating output directory"
    )
    
    # Add run target
    add_custom_target(run_poisson
        COMMAND ${CMAKE_BINARY_DIR}/bin/crossflow_poisson
        DEPENDS crossflow_poisson
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
endif()