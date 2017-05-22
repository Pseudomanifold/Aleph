# Optimizes for a useful debugging experience. While the test cases
# might run slightly faster if "real" optimizations were enabled, I
# like this flag for debugging.
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Og" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-g" )

#
# The tests should always have these nasty flags enabled. What's the
# use of the test without them?
#
# Coincidentally, isn't it nice how "-Wall" will not enable any of
# the flags described afterwards?
#
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wall" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wextra" )

# There are some issues with user-installed libraries in Mac OS X that
# cause warnings to appear in code that does *not* belong to Aleph, so
# I rather want this switch to be disabled.
#
# The automated tests under Linux will still have this setting though,
# so I should not miss to many issues.
IF( NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Werror" )
ENDIF()

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wconversion" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wnon-virtual-dtor" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wold-style-cast" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Woverloaded-virtual" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wself-init" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wunsafe-loop-optimization" )

IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
  IF( CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.4" )
    # Disable this warning for older versions of GCC because it will
    # be triggered by some of the test cases.
    SET_SOURCE_FILES_PROPERTIES(
      test_io_gml.cc
      PROPERTIES COMPILE_FLAGS "-Wno-missing-field-initializers"
    )
  ENDIF()
ENDIF()

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-pedantic" )

ADD_EXECUTABLE( test_boundary_matrix_reduction    test_boundary_matrix_reduction.cc )
ADD_EXECUTABLE( test_clique_enumeration           test_clique_enumeration.cc )
ADD_EXECUTABLE( test_clique_graph                 test_clique_graph.cc )
ADD_EXECUTABLE( test_connected_components         test_connected_components.cc )
ADD_EXECUTABLE( test_data_descriptors             test_data_descriptors.cc )
ADD_EXECUTABLE( test_filesystem                   test_filesystem.cc )
ADD_EXECUTABLE( test_graph_generation             test_graph_generation.cc )
ADD_EXECUTABLE( test_io_functions                 test_io_functions.cc )
ADD_EXECUTABLE( test_io_gml                       test_io_gml.cc )
ADD_EXECUTABLE( test_io_pajek                     test_io_pajek.cc )
ADD_EXECUTABLE( test_io_vtk                       test_io_vtk.cc )
ADD_EXECUTABLE( test_mesh                         test_mesh.cc )
ADD_EXECUTABLE( test_munkres                      test_munkres.cc )
ADD_EXECUTABLE( test_nearest_neighbours           test_nearest_neighbours.cc )
ADD_EXECUTABLE( test_persistence_diagrams         test_persistence_diagrams.cc )
ADD_EXECUTABLE( test_persistent_homology_complete test_persistent_homology_complete.cc )
ADD_EXECUTABLE( test_point_clouds                 test_point_clouds.cc )
ADD_EXECUTABLE( test_rips_expansion               test_rips_expansion.cc )
ADD_EXECUTABLE( test_rips_skeleton                test_rips_skeleton.cc )
ADD_EXECUTABLE( test_union_find                   test_union_find.cc )
ADD_EXECUTABLE( test_step_function                test_step_function.cc )
ADD_EXECUTABLE( test_wasserstein                  test_wasserstein.cc )

ADD_TEST( boundary_matrix_reduction    test_boundary_matrix_reduction )
ADD_TEST( clique_enumeration           test_clique_enumeration )
ADD_TEST( clique_graph                 test_clique_graph )
ADD_TEST( connected_components         test_connected_components )
ADD_TEST( data_descriptors             test_data_descriptors )
ADD_TEST( filesystem                   test_filesystem )
ADD_TEST( graph_generation             test_graph_generation )
ADD_TEST( io_functions                 test_io_functions )
ADD_TEST( io_gml                       test_io_gml )
ADD_TEST( io_pajek                     test_io_pajek )
ADD_TEST( io_vtk                       test_io_vtk )
ADD_TEST( mesh                         test_mesh )
ADD_TEST( munkres                      test_munkres )
ADD_TEST( nearest_neighbours           test_nearest_neighbours )
ADD_TEST( persistence_diagrams         test_persistence_diagrams )
ADD_TEST( persistent_homology_complete test_persistent_homology_complete )
ADD_TEST( point_clouds                 test_point_clouds )
ADD_TEST( rips_expansion               test_rips_expansion )
ADD_TEST( rips_skeleton                test_rips_skeleton )
ADD_TEST( step_function                test_step_function )
ADD_TEST( union_find                   test_union_find )
ADD_TEST( wasserstein                  test_wasserstein )
