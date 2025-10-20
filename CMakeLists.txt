cmake_minimum_required(VERSION 3.10)
project(LWproj)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")

# Find Boost Iostreams for numerical integration and OMP for parallelization
find_package(Boost REQUIRED COMPONENTS iostreams)
find_package(OpenMP REQUIRED)


# Add executable
add_executable(lwproj
    src/main.cpp
    src/plane_parallel.cpp
    src/MC.cpp
)

# Include directories
target_include_directories(lwproj PRIVATE 
    ${CMAKE_SOURCE_DIR}/include   # your gnuplot-iostream.h
    ${Boost_INCLUDE_DIRS}         # Boost headers
)

# Link Boost libraries
target_link_libraries(lwproj PRIVATE ${Boost_LIBRARIES})
target_link_libraries(lwproj PRIVATE OpenMP::OpenMP_CXX)
