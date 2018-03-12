function(setup_test test_file)
  file(GLOB source ${test_file}.c test.c)
	add_executable(${test_file} ${source})
	target_link_libraries(${test_file} mtmetis m)
	add_test(NAME ${test_file} COMMAND ${test_file})
endfunction()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/src/test)

include_directories(.)

setup_test(base_test)
setup_test(graph_test)
