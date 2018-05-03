IF( BUILD_TOOLS )
  MESSAGE( STATUS "Building tools" )

  ADD_EXECUTABLE( betti_numbers                                  betti_numbers.cc )
  ADD_EXECUTABLE( clique_communities_to_json                     clique_communities_to_json.cc )
  ADD_EXECUTABLE( clique_persistence_diagram                     clique_persistence_diagram.cc )
  ADD_EXECUTABLE( function_analysis                              function_analysis.cc )
  ADD_EXECUTABLE( graph_analysis                                 graph_analysis.cc )
  ADD_EXECUTABLE( interlevel_set_persistence_hierarchy           interlevel_set_persistence_hierarchy.cc )
  ADD_EXECUTABLE( local_dimensionality                           local_dimensionality.cc )
  ADD_EXECUTABLE( make_signature                                 make_signature.cc )
  ADD_EXECUTABLE( mean_curvature                                 mean_curvature.cc )
  ADD_EXECUTABLE( multi_scale_skeleton                           multi_scale_skeleton.cc )
  ADD_EXECUTABLE( persistence_diagram_entropies                  persistence_diagram_entropies.cc )
  ADD_EXECUTABLE( persistence_diagram_statistics                 persistence_diagram_statistics.cc )
  ADD_EXECUTABLE( persistence_indicator_function                 persistence_indicator_function.cc )
  ADD_EXECUTABLE( persistence_indicator_function_confidence_sets persistence_indicator_function_confidence_sets.cc )
  ADD_EXECUTABLE( persistence_indicator_function_glyph           persistence_indicator_function_glyph.cc )
  ADD_EXECUTABLE( persistent_intersection_homology               persistent_intersection_homology.cc )
  ADD_EXECUTABLE( persistent_spines                              persistent_spines.cc )
  ADD_EXECUTABLE( piecewise_linear_function_smoothing            piecewise_linear_function_smoothing.cc )
  ADD_EXECUTABLE( point_cloud_data_descriptors                   point_cloud_data_descriptors.cc )
  ADD_EXECUTABLE( process_image_patches                          process_image_patches.cc )
  ADD_EXECUTABLE( random_stratified_graphs                       random_stratified_graphs.cc )
  ADD_EXECUTABLE( relevant_edges                                 relevant_edges.cc )
  ADD_EXECUTABLE( rips_distributions                             rips_distributions.cc )
  ADD_EXECUTABLE( sparse_adjacency_matrices                      sparse_adjacency_matrices.cc )
  ADD_EXECUTABLE( spectral_analysis                              spectral_analysis.cc )
  ADD_EXECUTABLE( stratified_graph_analysis                      stratified_graph_analysis.cc )
  ADD_EXECUTABLE( topological_distance                           topological_distance.cc )
  ADD_EXECUTABLE( wicked_triangulations                          wicked_triangulations.cc )

  ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-O3" )

  ######################################################################
  # Check compiler flags
  ######################################################################
  #
  # A few targets require a modern compiler, so we need to make sure
  # that they are *not* built when an old compiler has been detected
  # or the flag for C++14 is not supported.

  CHECK_CXX_COMPILER_FLAG( "-std=c++14" ALEPH_HAVE_FLAG_CXX14 )

  IF( ALEPH_HAVE_FLAG_CXX14 )
    # The confidence set calculation requires the availability of generic
    # lambdas. Ditto for persistent intersection homology and smoothing.
    SET_SOURCE_FILES_PROPERTIES(
      persistence_indicator_function_confidence_sets.cc
      persistent_intersection_homology.cc
      piecewise_linear_function_smoothing.cc
      PROPERTIES COMPILE_FLAGS "-std=c++14"
    )

  # Old or unsupported compiler; let's exclude the targets from the
  # build process.
  ELSE()
    SET_TARGET_PROPERTIES(
      persistence_indicator_function_confidence_sets
      persistent_intersection_homology
      piecewise_linear_function_smoothing
      PROPERTIES EXCLUDE_FROM_ALL TRUE
    )
  ENDIF()

  IF( FLANN_FOUND AND FLANN_VERSION VERSION_GREATER "1.8.4" )
    TARGET_LINK_LIBRARIES( local_dimensionality             lz4 )
    TARGET_LINK_LIBRARIES( mean_curvature                   lz4 )
    TARGET_LINK_LIBRARIES( persistence_diagram_entropies    lz4 )
    TARGET_LINK_LIBRARIES( persistent_intersection_homology lz4 )
    TARGET_LINK_LIBRARIES( point_cloud_data_descriptors     lz4 )
    TARGET_LINK_LIBRARIES( rips_distributions               lz4 )
  ENDIF()

  IF( ALEPH_WITH_HDF5 )
    TARGET_INCLUDE_DIRECTORIES( relevant_edges
      SYSTEM PUBLIC
        ${HDF5_INCLUDE_DIRS}
        ${HDF5_CXX_INCLUDE_DIRS}
    )
    TARGET_LINK_LIBRARIES( relevant_edges
      ${HDF5_CXX_LIBRARIES}
      ${HDF5_LIBRARIES}
    )
  ENDIF()

ELSE()
  MESSAGE( STATUS "Not building tools (toggle BUILD_TOOLS to change this)" )
ENDIF()
