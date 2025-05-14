cmake_minimum_required(VERSION 3.15)

project(iSpan LANGUAGES CXX)
option(BUILD_MPI "Whether to build the iSpan MPI version" OFF)

# cxx standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# normalise output directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# enable lto when available
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
if (ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

# fetch required dependencies
find_package(OpenMP REQUIRED)
if (BUILD_MPI)
    find_package(MPI REQUIRED CXX)
endif ()

# add subdirectories
add_subdirectory(src)
if (BUILD_MPI)
    add_subdirectory(src_mpi)
endif ()
