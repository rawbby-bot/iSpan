cmake_minimum_required(VERSION 3.15)

add_executable(txt_to_bin_fw txt_to_bin_fw.cpp)
add_executable(txt_to_bin_bw txt_to_bin_bw.cpp)

add_custom_command(TARGET txt_to_bin_fw
        POST_BUILD
        COMMAND $<TARGET_FILE:txt_to_bin_fw> "${CMAKE_CURRENT_LIST_DIR}/toy_edgelist.txt"
        WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
        COMMENT "Test: running toy graph conversion post build")

add_custom_command(TARGET txt_to_bin_bw
        POST_BUILD
        COMMAND $<TARGET_FILE:txt_to_bin_bw> "${CMAKE_CURRENT_LIST_DIR}/toy_edgelist.txt"
        WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
        COMMENT "Test: running toy graph conversion post build")
