ADD_EXECUTABLE( clique_communities_to_json                     clique_communities_to_json.cc )
ADD_EXECUTABLE( clique_persistence_diagram                     clique_persistence_diagram.cc )
ADD_EXECUTABLE( interlevel_set_persistence_hierarchy           interlevel_set_persistence_hierarchy.cc )
ADD_EXECUTABLE( persistence_diagram_statistics                 persistence_diagram_statistics.cc )
ADD_EXECUTABLE( persistence_indicator_function                 persistence_indicator_function.cc )
ADD_EXECUTABLE( persistence_indicator_function_confidence_sets persistence_indicator_function_confidence_sets.cc )
ADD_EXECUTABLE( persistence_indicator_function_glyph           persistence_indicator_function_glyph.cc )
ADD_EXECUTABLE( point_cloud_data_descriptors                   point_cloud_data_descriptors.cc )
ADD_EXECUTABLE( sparse_adjacency_matrices                      sparse_adjacency_matrices.cc )
ADD_EXECUTABLE( topological_distance                           topological_distance.cc )
ADD_EXECUTABLE( wicked_triangulations                          wicked_triangulations.cc )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-O3" )

# The confidence set calculation requires the availability of generic
# lambdas.
SET_SOURCE_FILES_PROPERTIES(
  persistence_indicator_function_confidence_sets.cc
  PROPERTIES COMPILE_FLAGS "-std=c++14"
)

