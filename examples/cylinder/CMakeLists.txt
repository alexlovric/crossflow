# Add the cylinder example executable
add_executable(crossflow_cylinder
    main.cpp
)

# Link against the crossflow library and VTK
target_link_libraries(crossflow_cylinder
    PRIVATE
        crossflow
        ${VTK_LIBS}
)

# Set output directory
set_target_properties(crossflow_cylinder PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Create output directory and copy mesh file
add_custom_command(TARGET crossflow_cylinder POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin/out
    COMMAND ${CMAKE_COMMAND} -E copy_if_different 
        ${CMAKE_CURRENT_SOURCE_DIR}/cylinder.msh 
        ${CMAKE_BINARY_DIR}/bin/
    COMMENT "Setting up cylinder example files..."
    VERBATIM
)

# Add a custom target to run the example
add_custom_target(run_cylinder
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin/out
    COMMAND ${CMAKE_COMMAND} -E copy_if_different 
        ${CMAKE_CURRENT_SOURCE_DIR}/cylinder.msh 
        ${CMAKE_BINARY_DIR}/bin/
    COMMAND ./crossflow_cylinder
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    DEPENDS crossflow_cylinder
    COMMENT "Running cylinder example..."
)
