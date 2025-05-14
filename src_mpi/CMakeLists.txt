cmake_minimum_required(VERSION 3.15)

set(TARGET_NAME "${PROJECT_NAME}_mpi")

file(GLOB TARGET_SRC ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

add_executable(${TARGET_NAME} ${TARGET_SRC})

target_include_directories(${TARGET_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${TARGET_NAME} PRIVATE
        MPI::MPI_CXX
        OpenMP::OpenMP_CXX)
