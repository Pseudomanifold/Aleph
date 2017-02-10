# The tests should always have these nasty flags enabled. What's the
# use of the test without them?
#
# Coincidentally, isn't it nice how "-Wall" will not enable any of
# the flags described afterwards?
#
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wall" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wextra" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Werror" )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-WConversion" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wnon-virtual-dtor" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wold-style-cast" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Woverloaded-virtual" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wself-init" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wunsafe-loop-optimization" )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-pedantic" )

ADD_EXECUTABLE( test_boundary_matrix_reduction    test_boundary_matrix_reduction.cc )
ADD_EXECUTABLE( test_clique_graph                 test_clique_graph.cc )
ADD_EXECUTABLE( test_connected_components         test_connected_components.cc )
ADD_EXECUTABLE( test_data_descriptors             test_data_descriptors.cc )
ADD_EXECUTABLE( test_io_gml                       test_io_gml.cc )
ADD_EXECUTABLE( test_munkres                      test_munkres.cc )
ADD_EXECUTABLE( test_nearest_neighbours           test_nearest_neighbours.cc )
ADD_EXECUTABLE( test_persistent_homology_complete test_persistent_homology_complete.cc )
ADD_EXECUTABLE( test_point_clouds                 test_point_clouds.cc )
ADD_EXECUTABLE( test_rips_expansion               test_rips_expansion.cc )
ADD_EXECUTABLE( test_rips_skeleton                test_rips_skeleton.cc )
ADD_EXECUTABLE( test_step_function                test_step_function.cc )
ADD_EXECUTABLE( test_wasserstein                  test_wasserstein.cc )

ADD_TEST( boundary_matrix_reduction    test_boundary_matrix_reduction )
ADD_TEST( clique_graph                 test_clique_graph )
ADD_TEST( connected_components         test_connected_components )
ADD_TEST( data_descriptors             test_data_descriptors )
ADD_TEST( io_gml                       test_io_gml )
ADD_TEST( munkres                      test_munkres )
ADD_TEST( nearest_neighbours           test_nearest_neighbours )
ADD_TEST( persistent_homology_complete test_persistent_homology_complete )
ADD_TEST( point_clouds                 test_point_clouds )
ADD_TEST( rips_expansion               test_rips_expansion )
ADD_TEST( rips_skeleton                test_rips_skeleton )
ADD_TEST( step_function                test_step_function )
ADD_TEST( wasserstein                  test_wasserstein )
