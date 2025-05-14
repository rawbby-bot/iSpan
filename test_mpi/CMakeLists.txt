cmake_minimum_required(VERSION 3.15)

file(GLOB CPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach (TARGET_SRC IN LISTS CPP_FILES)
    get_filename_component(TARGET_NAME ${TARGET_SRC} NAME_WE)

    add_executable(${TARGET_NAME} ${TARGET_SRC})
    add_dependencies(${TARGET_NAME} txt_to_bin_fw)
    add_dependencies(${TARGET_NAME} txt_to_bin_bw)

    target_link_libraries(${TARGET_NAME} PRIVATE
            MPI::MPI_CXX
            OpenMP::OpenMP_CXX)

    add_test(NAME ${TARGET_NAME}
            COMMAND mpirun -np 2 $<TARGET_FILE:${TARGET_NAME}>
            WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

endforeach ()

# as src and src_mpi are not libraries but monolithic CLIs we 'steal' implementations rather than linking
set(EXTRA_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/../src_mpi/graph.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/../src_mpi/load_print.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/../src_mpi/scc_core.cpp")
target_sources(test_toy_mpi PRIVATE ${EXTRA_SOURCES})
target_sources(test_toy2_mpi PRIVATE ${EXTRA_SOURCES})
add_dependencies(test_toy_mpi convert_toy)
add_dependencies(test_toy2_mpi convert_toy2)
