# Define the library and its source files
add_library(NSE_core_lib STATIC
                
)

# Tell CMake where the headers (like LGFCore.H) are
target_include_directories(NSE_core_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link to AMReX so the core knows about MultiFabs and MPI
target_link_libraries(NSE_core_lib PUBLIC AMReX::amrex)