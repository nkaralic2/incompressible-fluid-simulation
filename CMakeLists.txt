cmake_minimum_required(VERSION 3.16)
project(FluidSimulation LANGUAGES CXX)

enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release by default for benchmarks
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler optimization flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
elif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /W4")
endif()

# Include directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files for the main simulator library
set(SIM_CORE_SOURCES
    src/Grid.cpp
    src/FLIPSimulator.cpp
    src/PressureSolver.cpp
)

set(SIM_SOURCES
    ${SIM_CORE_SOURCES}
    src/Benchmark.cpp
)

# Main executable
add_executable(fluid_sim
    src/main.cpp
    ${SIM_SOURCES}
)

add_test(NAME fluid_sim_correctness COMMAND fluid_sim --test)

option(BUILD_SFML_VIEWER "Build the real-time SFML visualization" ON)
if(BUILD_SFML_VIEWER)
    find_package(SFML 3 COMPONENTS Graphics Window System CONFIG QUIET)
    if(SFML_FOUND)
        add_executable(fluid_viewer
            src/SFMLViewer.cpp
            ${SIM_CORE_SOURCES}
        )
        target_compile_definitions(fluid_viewer PRIVATE
            WATER_SHADER_PATH="${CMAKE_SOURCE_DIR}/shaders/water.frag"
        )
        target_link_libraries(fluid_viewer PRIVATE SFML::Graphics SFML::Window SFML::System)
        message(STATUS "SFML viewer: enabled")
    else()
        message(WARNING "SFML not found: fluid_viewer will not be built")
    endif()
endif()

# Create output directories
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/benchmarks)

# Print configuration
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ standard: C++${CMAKE_CXX_STANDARD}")
