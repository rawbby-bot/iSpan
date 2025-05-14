cmake_minimum_required(VERSION 3.15)

add_executable(txt_to_bin_fw txt_to_bin_fw.cpp)
add_executable(txt_to_bin_bw txt_to_bin_bw.cpp)

function(convert_graph_post_build BASE_DIR CODE)

    set(TEMP_DIR "${CMAKE_BINARY_DIR}/temp/${CODE}")
    file(MAKE_DIRECTORY "${TEMP_DIR}")

    set(OUTPUT_FILES
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_fw_adjacent.bin"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_fw_begin.bin"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_fw_head.bin"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_out_degree.bin"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_bw_adjacent.bin"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_bw_begin.bin"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_bw_head.bin"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_in_degree.bin")

    add_custom_command(
            OUTPUT ${OUTPUT_FILES}
            COMMAND $<TARGET_FILE:txt_to_bin_fw> "${BASE_DIR}/${CODE}_edgelist.txt"
            COMMAND $<TARGET_FILE:txt_to_bin_bw> "${BASE_DIR}/${CODE}_edgelist.txt"
            COMMAND ${CMAKE_COMMAND} -E rename "fw_adjacent.bin" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_fw_adjacent.bin"
            COMMAND ${CMAKE_COMMAND} -E rename "fw_begin.bin" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_fw_begin.bin"
            COMMAND ${CMAKE_COMMAND} -E rename "fw_head.bin" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_fw_head.bin"
            COMMAND ${CMAKE_COMMAND} -E rename "out_degree.bin" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_out_degree.bin"
            COMMAND ${CMAKE_COMMAND} -E rename "bw_adjacent.bin" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_bw_adjacent.bin"
            COMMAND ${CMAKE_COMMAND} -E rename "bw_begin.bin" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_bw_begin.bin"
            COMMAND ${CMAKE_COMMAND} -E rename "bw_head.bin" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_bw_head.bin"
            COMMAND ${CMAKE_COMMAND} -E rename "in_degree.bin" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CODE}_in_degree.bin"
            DEPENDS txt_to_bin_bw txt_to_bin_fw
            WORKING_DIRECTORY "${TEMP_DIR}")

    add_custom_target("convert_${CODE}" ALL DEPENDS ${OUTPUT_FILES})

endfunction()

convert_graph_post_build("${CMAKE_CURRENT_LIST_DIR}" "toy")
convert_graph_post_build("${CMAKE_CURRENT_LIST_DIR}" "toy2")
