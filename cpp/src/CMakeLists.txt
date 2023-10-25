find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(MPI REQUIRED)
find_package(MAGMA REQUIRED)

add_executable(hf hf.cpp input.cpp helper.cpp)
target_include_directories(hf PUBLIC ${CMAKE_SOURCE_DIR}/include
    /usr/local/include
    /usr/local/include/eigen3
     /usr/include/eigen3
     /usr/include
    )
