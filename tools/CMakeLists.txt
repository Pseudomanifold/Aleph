ADD_EXECUTABLE( clique-dna                     clique-dna.cc )
ADD_EXECUTABLE( clique-persistence-diagram     clique-persistence-diagram.cc )
ADD_EXECUTABLE( persistence-indicator-function persistence-indicator-function.cc )
ADD_EXECUTABLE( wasserstein-distance           wasserstein-distance.cc )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-O3" )
