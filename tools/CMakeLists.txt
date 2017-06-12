ADD_EXECUTABLE( create_functions                        create_functions.cc )
ADD_EXECUTABLE( clique-communities                      clique-communities.cc )
ADD_EXECUTABLE( clique-dna                              clique-dna.cc )
ADD_EXECUTABLE( clique-persistence-diagram              clique-persistence-diagram.cc )
ADD_EXECUTABLE( extended_persistence_hierarchy          extended_persistence_hierarchy.cc )
ADD_EXECUTABLE( network-persistence-diagrams            network-persistence-diagrams.cc )

ADD_EXECUTABLE( persistence_diagram_statistics          persistence_diagram_statistics.cc )
ADD_EXECUTABLE( persistence_indicator_function          persistence_indicator_function.cc )
ADD_EXECUTABLE( topological_distance                    topological_distance.cc )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-O3" )
