# magma-config.cmake

# Get the Conda environment path from the environment variable
set(CONDA_ENV_PATH $ENV{CONDA_PREFIX})

# Set up the MAGMA variables pointing to the Conda directories
set(MAGMA_INCLUDE_DIR "${CONDA_ENV_PATH}/include")
set(MAGMA_LIBRARY "${CONDA_ENV_PATH}/lib/libmagma.so")  # Modify this as necessary

# Check if the directories and files exist
if(NOT EXISTS "${MAGMA_INCLUDE_DIR}")
  message(FATAL_ERROR "Could not find MAGMA include directory")
endif()
if(NOT EXISTS "${MAGMA_LIBRARY}")
  message(FATAL_ERROR "Could not find MAGMA library")
endif()

# Attempt to find CUDA
find_package(CUDA QUIET)

# Create the MAGMA::MAGMA target (could be CPU-only or with CUDA)
if(NOT TARGET MAGMA::MAGMA)
  add_library(MAGMA::MAGMA INTERFACE IMPORTED)
endif()

if(CUDA_FOUND)
  message(STATUS "Found CUDA: ${CUDA_INCLUDE_DIRS}")
  set_target_properties(MAGMA::MAGMA PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MAGMA_INCLUDE_DIR};${CUDA_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${MAGMA_LIBRARY};${CUDA_LIBRARIES}"
  )
else()
message(STATUS "Could not find CUDA, using CPU-only MAGMA")
  set_target_properties(MAGMA::MAGMA PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MAGMA_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${MAGMA_LIBRARY}"
  )
endif()

