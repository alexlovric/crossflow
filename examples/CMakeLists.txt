cmake_minimum_required(VERSION 3.10)
project(CrossFlowExamples VERSION 1.0.0 LANGUAGES CXX)

# Standardize on C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Link against main library and VTK
# add_subdirectory(poisson)
add_subdirectory(cylinder)
