# Where the header files reside
#include_directories(../)

# Build program.
add_executable(strings strings.c)
add_executable(gksort gksort.c)
add_executable(fis fis.c)
add_executable(rw rw.c)
add_executable(gkgraph gkgraph.c)
add_executable(csrcnv csrcnv.c)
add_executable(grKx grKx.c)
foreach(prog strings gksort fis rw gkgraph csrcnv grKx)
  target_link_libraries(${prog} GKlib)
endforeach(prog)

