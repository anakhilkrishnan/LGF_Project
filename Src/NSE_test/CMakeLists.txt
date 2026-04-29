add_executable(nserun
                main.cpp)

# uncomment below line for gpu compilation
setup_target_for_cuda_compilation(nserun)

# Link to the SAME LGF_core_lib to use MeshData and kernels
target_link_libraries(nserun PRIVATE LGF_core_lib NSE_core_lib AMReX::amrex)

set_target_properties(nserun PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/Exec/NSE_test"
)

# Define the path to your inputs file (assumes it's in the same folder as CMakeLists)
set(INPUTS_FILE "${CMAKE_CURRENT_SOURCE_DIR}/inputs")

# Update the input file sync command to follow the binary
add_custom_command(TARGET nserun POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_SOURCE_DIR}/inputs"
    "${CMAKE_SOURCE_DIR}/Exec/NSE_test/inputs"
    COMMENT "Syncing inputs to Exec/NSE_test..."
)

# Update the jobscript file sync command to follow the binary
add_custom_command(TARGET nserun POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_SOURCE_DIR}/debug.sh"
    "${CMAKE_SOURCE_DIR}/Exec/NSE_test/debug.sh"
    COMMENT "Syncing jobscript to Exec/NSE_test..."
)