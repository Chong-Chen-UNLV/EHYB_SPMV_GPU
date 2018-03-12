function(setup_test test_file)
  file(GLOB source ${test_file}.cpp)
	add_executable(${test_file} ${source})
	target_link_libraries(${test_file} wildriver)
	add_test(NAME ${test_file} COMMAND ${test_file})
endfunction()


include_directories(.)


setup_test(CSRFile_test)
setup_test(GraphInHandle_test)
setup_test(GraphOutHandle_test)
setup_test(MatrixInHandle_test)
setup_test(MatrixOutHandle_test)
setup_test(MetisFile_test)
setup_test(WildRiver_test)
setup_test(PlainVectorFile_test)

