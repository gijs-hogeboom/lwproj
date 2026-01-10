cmake_minimum_required(VERSION 3.10)
project(LWproj)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find OpenMP for parallelization
find_package(OpenMP REQUIRED)

# NetCDF include path
execute_process(
    COMMAND nc-config --cflags
    OUTPUT_VARIABLE NETCDF_CFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# NetCDF libraries
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

# Include directories
target_include_directories(lwproj PRIVATE 
    ${CMAKE_SOURCE_DIR}/include   # e.g. gnuplot-iostream.h
)

# Add NetCDF compile options
target_compile_options(lwproj PRIVATE ${NETCDF_CFLAGS})

# Link OpenMP and NetCDF libraries
target_link_libraries(lwproj PRIVATE OpenMP::OpenMP_CXX ${NETCDF_LIBS} netcdf_c++4)

# ---------------------------
# Debug build: enable sanitizers
# ---------------------------
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Debug build: enabling AddressSanitizer")

    # AddressSanitizer flags (safe with GCC + OpenMP)
    target_compile_options(lwproj PRIVATE
        -g
        -O1
        -fsanitize=address
        -fno-omit-frame-pointer
        -fno-inline
        -fopenmp
    )

    target_link_options(lwproj PRIVATE
        -fsanitize=address
        -fopenmp
    )

    # If you want ThreadSanitizer instead (Clang + OpenMP), comment the above and use:
    # target_compile_options(lwproj PRIVATE
    #     -g -O1 -fsanitize=thread -fno-omit-frame-pointer -fno-inline -fopenmp
    # )
    # target_link_options(lwproj PRIVATE -fsanitize=thread -fopenmp)
endif()

