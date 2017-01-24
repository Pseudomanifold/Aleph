ADD_EXECUTABLE( wasserstein-distance wasserstein-distance.cc )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-O3" )
