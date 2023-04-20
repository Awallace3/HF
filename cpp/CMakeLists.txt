cmake_minimum_required(VERSION 3.10)

project(HF VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++-12")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fopenmp -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -fopenmp")

include_directories(/usr/local/include)
link_directories(/usr/local/lib)
# target_include_directories(hf PUBLIC include)

# add_library(eigen3 INTERFACE)
# add_library(iostream INTERFACE)

find_package(OpenMP REQUIRED)
# find_package(MPI REQUIRED)

# add_subdirectory(include)
add_subdirectory(src)

# add_executable(hf main.cpp)
