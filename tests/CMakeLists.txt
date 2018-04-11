CONFIGURE_FILE( Base.hh.in ${CMAKE_SOURCE_DIR}/tests/Base.hh )

IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
  IF( CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.4.1" )
    # Disable this warning for older versions of GCC because it will
    # be triggered by some of the test cases.
    SET_SOURCE_FILES_PROPERTIES(
      test_io_pajek.cc
      test_io_gml.cc
      test_io_vtk.cc
      PROPERTIES COMPILE_FLAGS "-Wno-missing-field-initializers"
    )

    # Disable warnings for unknown pragmas in case OpenMP cannot be
    # found in the specified configuration. This will only make the
    # test take longer.
    SET_SOURCE_FILES_PROPERTIES(
      test_persistence_diagrams.cc
      PROPERTIES COMPILE_FLAGS "-Wno-unknown-pragmas"
    )
  ENDIF()
ENDIF()

IF( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  # Ignore warnings related to the GNU statement expressions if we
  # compile with clang.
  ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wno-gnu-statement-expression" )
ENDIF()

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
# so I should not miss too many issues.
IF( NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Werror" )
ENDIF()

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wconversion" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wnon-virtual-dtor" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wold-style-cast" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Woverloaded-virtual" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wself-init" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wunsafe-loop-optimization" )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-pedantic" )

# FIXME: Ugly hack to make generic lambdas work
SET_SOURCE_FILES_PROPERTIES(
  test_bootstrap.cc
  PROPERTIES COMPILE_FLAGS "-std=c++14"
)

ADD_EXECUTABLE( test_barycentric_subdivision          test_barycentric_subdivision.cc )
ADD_EXECUTABLE( test_beta_skeleton                    test_beta_skeleton.cc )
ADD_EXECUTABLE( test_bootstrap                        test_bootstrap.cc )
ADD_EXECUTABLE( test_boundary_matrix_reduction        test_boundary_matrix_reduction.cc )
ADD_EXECUTABLE( test_cech_expansion                   test_cech_expansion.cc )
ADD_EXECUTABLE( test_clique_enumeration               test_clique_enumeration.cc )
ADD_EXECUTABLE( test_clique_graph                     test_clique_graph.cc )
ADD_EXECUTABLE( test_combinatorial_curvature          test_combinatorial_curvature.cc )
ADD_EXECUTABLE( test_connected_components             test_connected_components.cc )
ADD_EXECUTABLE( test_cover_tree                       test_cover_tree.cc )
ADD_EXECUTABLE( test_data_descriptors                 test_data_descriptors.cc )
ADD_EXECUTABLE( test_distances                        test_distances.cc )
ADD_EXECUTABLE( test_dowker_complex                   test_dowker_complex.cc )
ADD_EXECUTABLE( test_filesystem                       test_filesystem.cc )
ADD_EXECUTABLE( test_fractal_dimension                test_fractal_dimension.cc )
ADD_EXECUTABLE( test_graph_generation                 test_graph_generation.cc )
ADD_EXECUTABLE( test_floyd_warshall                   test_floyd_warshall.cc )
ADD_EXECUTABLE( test_heat_kernel                      test_heat_kernel.cc )
ADD_EXECUTABLE( test_io_bipartite_adjacency_matrix    test_io_bipartite_adjacency_matrix.cc )
ADD_EXECUTABLE( test_io_functions                     test_io_functions.cc )
ADD_EXECUTABLE( test_io_gml                           test_io_gml.cc )
ADD_EXECUTABLE( test_io_graphml                       test_io_graphml.cc )
ADD_EXECUTABLE( test_io_hdf5                          test_io_hdf5.cc )
ADD_EXECUTABLE( test_io_json                          test_io_json.cc )
ADD_EXECUTABLE( test_io_lexicographic_triangulation   test_io_lexicographic_triangulation.cc )
ADD_EXECUTABLE( test_io_pajek                         test_io_pajek.cc )
ADD_EXECUTABLE( test_io_sparse_adjacency_matrix       test_io_sparse_adjacency_matrix.cc )
ADD_EXECUTABLE( test_io_vtk                           test_io_vtk.cc )
ADD_EXECUTABLE( test_kernel_density_estimator         test_kernel_density_estimator.cc )
ADD_EXECUTABLE( test_mesh                             test_mesh.cc )
ADD_EXECUTABLE( test_munkres                          test_munkres.cc )
ADD_EXECUTABLE( test_nearest_neighbours               test_nearest_neighbours.cc )
ADD_EXECUTABLE( test_partitions                       test_partitions.cc )
ADD_EXECUTABLE( test_persistence_diagrams             test_persistence_diagrams.cc )
ADD_EXECUTABLE( test_persistent_homology_complete     test_persistent_homology_complete.cc )
ADD_EXECUTABLE( test_persistent_intersection_homology test_persistent_intersection_homology.cc )
ADD_EXECUTABLE( test_piecewise_linear_function        test_piecewise_linear_function.cc )
ADD_EXECUTABLE( test_principal_component_analysis     test_principal_component_analysis.cc )
ADD_EXECUTABLE( test_point_clouds                     test_point_clouds.cc )
ADD_EXECUTABLE( test_rips_expansion                   test_rips_expansion.cc )
ADD_EXECUTABLE( test_rips_skeleton                    test_rips_skeleton.cc )
ADD_EXECUTABLE( test_spine                            test_spine.cc )
ADD_EXECUTABLE( test_tangent_space                    test_tangent_space.cc )
ADD_EXECUTABLE( test_union_find                       test_union_find.cc )
ADD_EXECUTABLE( test_step_function                    test_step_function.cc )
ADD_EXECUTABLE( test_witness_complex                  test_witness_complex.cc )

ADD_TEST( barycentric_subdivision          test_barycentric_subdivision )
ADD_TEST( beta_skeleton                    test_beta_skeleton )
ADD_TEST( bootstrap                        test_bootstrap )
ADD_TEST( boundary_matrix_reduction        test_boundary_matrix_reduction )
ADD_TEST( cech_expansion                   test_cech_expansion )
ADD_TEST( clique_enumeration               test_clique_enumeration )
ADD_TEST( clique_graph                     test_clique_graph )
ADD_TEST( combinatorial_curvature          test_combinatorial_curvature )
ADD_TEST( connected_components             test_connected_components )
ADD_TEST( data_descriptors                 test_data_descriptors )
ADD_TEST( distances                        test_distances )
ADD_TEST( dowker_complex                   test_dowker_complex )
ADD_TEST( filesystem                       test_filesystem )
ADD_TEST( fractal_dimension                test_fractal_dimension )
ADD_TEST( graph_generation                 test_graph_generation )
ADD_TEST( heat_kernel                      test_heat_kernel )
ADD_TEST( io_bipartite_adjacency_matrix    test_io_bipartite_adjacency_matrix )
ADD_TEST( io_functions                     test_io_functions )
ADD_TEST( io_gml                           test_io_gml )

# The test will build nonetheless, but the results will of course be
# incorrect if the library is not available.
IF( ALEPH_WITH_HDF5 )
  ADD_TEST( io_hdf5                        test_io_hdf5 )
ENDIF()

# The test will build nonetheless, but the results will of course be
# incorrect if the library is not available.
IF( ALEPH_WITH_RAPID_JSON )
  ADD_TEST( io_json                        test_io_json )
ENDIF()

ADD_TEST( io_lexicographic_triangulation   test_io_lexicographic_triangulation )
ADD_TEST( io_pajek                         test_io_pajek )
ADD_TEST( io_sparse_adjacency_matrix       test_io_sparse_adjacency_matrix )
ADD_TEST( io_vtk                           test_io_vtk )
ADD_TEST( kernel_density_estimator         test_kernel_density_estimator )
ADD_TEST( mesh                             test_mesh )
ADD_TEST( munkres                          test_munkres )
ADD_TEST( nearest_neighbours               test_nearest_neighbours )
ADD_TEST( partitions                       test_partitions )
ADD_TEST( persistence_diagrams             test_persistence_diagrams )
ADD_TEST( persistent_homology_complete     test_persistent_homology_complete )
ADD_TEST( persistent_intersection_homology test_persistent_intersection_homology )
ADD_TEST( piecewise_linear_function        test_piecewise_linear_function )
ADD_TEST( principal_component_analysis     test_principal_component_analysis )
ADD_TEST( point_clouds                     test_point_clouds )
ADD_TEST( rips_expansion                   test_rips_expansion )
ADD_TEST( rips_skeleton                    test_rips_skeleton )
ADD_TEST( spine                            test_spine )
ADD_TEST( step_function                    test_step_function )
ADD_TEST( tangent_space                    test_tangent_space )
ADD_TEST( union_find                       test_union_find )
ADD_TEST( witness_complex                  test_witness_complex )

# These test are a little bit special because they depend on another
# shared library so we need to adjust the interface if we are unable
# to build them.
IF( ALEPH_WITH_HDF5 )

  FOREACH( TARGET_NAME
    IN ITEMS
      test_io_gml
      test_io_graphml
      test_io_hdf5
      test_io_pajek
      test_io_vtk
  )
    TARGET_INCLUDE_DIRECTORIES( ${TARGET_NAME}
      SYSTEM PUBLIC
        ${HDF5_INCLUDE_DIRS}
        ${HDF5_CXX_INCLUDE_DIRS}
    )
    TARGET_LINK_LIBRARIES( ${TARGET_NAME}
      ${HDF5_CXX_LIBRARIES}
    )
  ENDFOREACH()
ENDIF()

IF( ALEPH_WITH_TINYXML2 )
  INCLUDE( ${tinyxml2_CONFIG} )
  TARGET_LINK_LIBRARIES( test_io_graphml tinyxml2 )
ENDIF()

IF( FLANN_FOUND AND FLANN_VERSION VERSION_GREATER "1.8.4" )
  TARGET_LINK_LIBRARIES( test_nearest_neighbours lz4 )
  TARGET_LINK_LIBRARIES( test_rips_skeleton      lz4 )
  TARGET_LINK_LIBRARIES( test_tangent_space      lz4 )
ENDIF()

########################################################################
# Python integration
########################################################################
#
# Special handling for the Python integration test. This program needs
# to be run but with an updated `PYTHONPATH`.

IF( BUILD_PYTHON_BINDINGS AND Python3_FOUND )
  # TODO: handle this dependency better by restructuring the way
  # bindings are being generated
  FIND_PACKAGE( pybind11 )

  IF( PYBIND11_FOUND )
    ADD_TEST(
      NAME python_integration
      COMMAND python3 test_python_integration.py
    )

    SET_TESTS_PROPERTIES( python_integration
      PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/bindings/python"
    )

    CONFIGURE_FILE( test_python_integration.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  ENDIF()
ENDIF()
