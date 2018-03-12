include_directories(.)
file(GLOB mtmetis_sources *.c)
file(GLOB domlib_sources 
  ${CMAKE_SOURCE_DIR}/${DOMLIB_PATH}/*.c)
file(GLOB metis_sources 
  ${METIS_LIB_PATH}/*.c
  ${METIS_LIB_PATH}/../GKlib/memory.c
  ${METIS_LIB_PATH}/../GKlib/mcore.c
  ${METIS_LIB_PATH}/../GKlib/timers.c
  ${METIS_LIB_PATH}/../GKlib/random.c
  ${METIS_LIB_PATH}/../GKlib/fs.c
  ${METIS_LIB_PATH}/../GKlib/error.c
  ${METIS_LIB_PATH}/../GKlib/io.c
  ${METIS_LIB_PATH}/../GKlib/gk_util.c
  ${METIS_LIB_PATH}/../GKlib/string.c
)

# build library
add_library(mtmetis ${MTMETIS_LIBRARY_TYPE} 
  ${mtmetis_sources}
  ${metis_sources}
  ${domlib_sources}
) 

# where to install libraries and binaries
install(TARGETS mtmetis 
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)

# binary
add_executable(mtmetis_bin mtmetis_bin.c)
set_target_properties(mtmetis_bin PROPERTIES OUTPUT_NAME mtmetis)
add_dependencies(mtmetis_bin wildriver)
add_dependencies(mtmetis_bin mtmetis)
target_link_libraries(mtmetis_bin mtmetis wildriver ${REQ_LIBS})
install(TARGETS mtmetis_bin
  RUNTIME DESTINATION bin
)


if (DEFINED TESTS AND NOT TESTS EQUAL 0)
  add_subdirectory("test")
endif()
