cmake_minimum_required(VERSION 3.10)
project(LWproj)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3 -fno-inline -fno-omit-frame-pointer")
set(CMAKE_BUILD_TYPE Release)

# Find Boost Iostreams for numerical integration and OMP for parallelization
find_package(Boost REQUIRED COMPONENTS iostreams)
find_package(OpenMP REQUIRED)

# NetCDF include path
execute_process(
	COMMAND nc-config --cflags
	OUTPUT_VARIABLE NETCDF_CFLAGS
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Netcdf libraries
execute_process(
	COMMAND nc-config --libs
	OUTPUT_VARIABLE NETCDF_LIBS
	OUTPUT_STRIP_TRAILING_WHITESPACE
)


# Add executable
add_executable(lwproj
    src/main.cpp
    src/plane_parallel.cpp
    src/MC.cpp
    src/photprop.cpp
)

target_compile_options(lwproj PRIVATE ${NETCDF_CFLAGS})

# Include directories
target_include_directories(lwproj PRIVATE 
    ${CMAKE_SOURCE_DIR}/include   # your gnuplot-iostream.h
    ${Boost_INCLUDE_DIRS}         # Boost headers
)

# Link Boost, OpenMP and NetCDF libraries
target_link_libraries(lwproj PRIVATE ${Boost_LIBRARIES})
target_link_libraries(lwproj PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(lwproj PRIVATE ${NETCDF_LIBS} netcdf_c++4)
