ADD_EXECUTABLE( network_analysis network_analysis.cc )
ADD_EXECUTABLE( ply              ply.cc              )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-O3" )
