# The tests should always have these nasty flags enabled. What's the
# use of the test without them?
#
# Coincidentally, isn't it nice how "-Wall" will not enable any of
# the flags described afterwards?
#
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wall" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wextra" )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-WConversion" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wnon-virtual-dtor" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wold-style-cast" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Woverloaded-virtual" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wself-init" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wunsafe-loop-optimization" )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-pedantic" )

ADD_EXECUTABLE( test_munkres            test_munkres.cc )
ADD_EXECUTABLE( test_nearest_neighbours test_nearest_neighbours.cc )
ADD_EXECUTABLE( test_point_clouds       test_point_clouds.cc )
ADD_EXECUTABLE( test_wasserstein        test_wasserstein.cc )

ADD_TEST( point_clouds test_point_clouds )
