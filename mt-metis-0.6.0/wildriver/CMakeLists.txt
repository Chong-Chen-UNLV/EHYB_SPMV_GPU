cmake_minimum_required(VERSION 2.8)
project(WILDRIVER)

if (DEFINED DEBUG)
  add_definitions(-DDEBUG=${DEBUG})
  set(CMAKE_BUILD_TYPE Debug)
  message("Building a debug release")
else()
  set(CMAKE_BUILD_TYPE Release)
  add_definitions(-DNDEBUG=1)
  message("Building an optimized release")
endif()

if (SHARED)
  set(WILDRIVER_LIBRARY_TYPE SHARED)
else()
  set(WILDRIVER_LIBRARY_TYPE STATIC)
endif()

if (DEFINED DEVEL AND NOT DEVEL EQUAL 0)
  message("Development warnings enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmissing-declarations") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-aliasing") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic") 
endif()

# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
			message("WARNING: Attempting to substitute -std=c++0x for -std=c++11.")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Doverride=") 
    else()
			# set c++11 for g++ >= 4.7
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 
    endif() 
else()
  # set c++11 standard for other compilers
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 
endif()


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# must be before including src
if (DEFINED TESTS AND NOT TESTS EQUAL 0)
  enable_testing()
endif()

include_directories("include")
add_subdirectory("include")
add_subdirectory("src")



