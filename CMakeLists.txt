
cmake_minimum_required(VERSION 3.16)

project(FluidSimulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(FluidSimulation
    source/main.cpp
    source/Grid.cpp
    source/FluidSimulation.cpp
)

target_include_directories(FluidSimulation PRIVATE source)
