include_directories(.)
file(GLOB sources *.cpp)

# libraries
add_library(wildriver ${WILDRIVER_LIBRARY_TYPE} 
  ${sources}
) 

install(TARGETS wildriver
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

if (DEFINED TESTS AND NOT TESTS EQUAL 0)
  add_subdirectory("test")
endif()
