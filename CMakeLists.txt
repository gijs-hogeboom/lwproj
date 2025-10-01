cmake_minimum_required(VERSION 3.10)
project(LWproj)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(include)

# Add executable and source files
add_executable(lwproj
    src/main.cpp
    src/plane_parallel.cpp
    src/MC.cpp
)
    
