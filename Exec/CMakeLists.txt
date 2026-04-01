cmake_minimum_required(VERSION 3.16)

project(LGF_Project
        LANGUAGES C CXX CUDA)

if(NOT DEFINED AMReX_ROOT)
  message("-- Download and configure AMReX from GitHub")

  #Download AMReX from GitHub
  include(FetchContent)
  set(FETCHCONTENT_QUIET OFF)

  FetchContent_Declare(
    amrex_code
    GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git/
    GIT_TAG        origin/development
    )

  FetchContent_MakeAvailable(amrex_code)


else()

  # Add AMReX
  message("-- Searching for AMReX install directory at ${AMReX_ROOT}")
  find_package(AMReX REQUIRED PATHS ${AMReX_ROOT}/lib/cmake/AMReX/AMReXConfig.cmake)

endif()

add_executable(lgfpoisson
                ../Src/main.cpp
                ../Src/createMultiFab.cpp
                ../Src/initializeMultiFab.cpp
                ../Src/addEverySourceBox.cpp
                ../Src/consolidateMultiFab.cpp
                ../Src/MyFunctions.H)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})



#setup_target_for_cuda_compilation(lgfpoisson)

target_include_directories(lgfpoisson PRIVATE ../Src)
target_link_libraries(lgfpoisson PRIVATE AMReX::amrex)

target_compile_definitions(lgfpoisson PRIVATE REL_RESULTS_DIR="Results/")

## Copy input files
file( COPY ${CMAKE_SOURCE_DIR}/inputs DESTINATION . )