cmake_minimum_required(VERSION 3.10)

project(HF VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fopenmp -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -fopenmp")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CONDA_ENV_PATH $ENV{CONDA_PREFIX})
message(STATUS "CONDA_ENV_PATH is set to ${CONDA_ENV_PATH}")
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR})

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(MPI REQUIRED)
find_package(MAGMA REQUIRED)
# find_package(CUDA REQUIRED)

# Check if the header and library were found
if(NOT MAGMA_INCLUDE_DIR)
  message(FATAL_ERROR "Could not find MAGMA include directory")
endif()
if(NOT MAGMA_LIBRARY)
  message(FATAL_ERROR "Could not find MAGMA library")
endif()

add_subdirectory(src)
