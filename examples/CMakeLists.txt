ADD_EXECUTABLE( 1D               1D.cc               )
ADD_EXECUTABLE( network_analysis network_analysis.cc )
ADD_EXECUTABLE( ply              ply.cc              )
ADD_EXECUTABLE( vietoris-rips    vietoris-rips.cc    )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-O3" )
