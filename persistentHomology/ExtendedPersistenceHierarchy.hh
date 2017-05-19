#ifndef ALEPH_PERSISTENT_HOMOLOGY_EXTENDED_PERSISTENCE_HIERARCHY__
#define ALEPH_PERSISTENT_HOMOLOGY_EXTENDED_PERSISTENCE_HIERARCHY__

#include <algorithm>
#include <map>
#include <stdexcept>
#include <set>
#include <vector>

#include <boost/bimap.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include "persistentHomology/PersistencePairing.hh"

#include "topology/SimplicialComplex.hh"
#include "topology/UnionFind.hh"

// FIXME: remove after debugging
#include <iostream>

namespace aleph
{

namespace detail
{

using AdjacencyGraph = boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::undirectedS,
  boost::no_property,
  boost::no_property>;

using SizeType = std::size_t;

/** Extracts the adjacency graph of a simplicial complex */
template <class Simplex>
std::pair<boost::bimap<typename Simplex::VertexType, SizeType>, AdjacencyGraph> extractZeroDimensionalAdjacencyGraph( const topology::SimplicialComplex<Simplex>& S )
{
  AdjacencyGraph adjacencyGraph;

  using Vertex           = typename Simplex::VertexType;
  using VertexDescriptor = boost::graph_traits<AdjacencyGraph>::vertex_descriptor;
  using bimap_type       = boost::bimap<Vertex, SizeType>;

  bimap_type indices;
  std::map<Vertex, VertexDescriptor> vdm;

  SizeType vertexIndex = 0;

  for( auto it = S.begin_dimension(); it != S.end_dimension(); ++it )
  {
    // Vertices
    if( it->dimension() == 0 )
    {
      indices.insert( typename bimap_type::value_type( *it->begin(), vertexIndex++ ) );

      auto vertex         = boost::add_vertex( adjacencyGraph );
      vdm[ *it->begin() ] = vertex;
    }

    // Edges
    else if( it->dimension() == 1 )
    {
      auto&& edge = *it;

      auto u = *( edge.begin()     );
      auto v = *( edge.begin() + 1 );

      // It is possible that the simplicial complex, being a _part_ of
      // a larger filtration, contains edges for which no vertices are
      // available.
      if( S.contains(u) && S.contains(v) )
      {
        boost::add_edge( vdm.at(u),
                         vdm.at(v),
                         adjacencyGraph );
      }
    }
  }

  return std::make_pair( indices, adjacencyGraph );
}

} // namespace detail

/** Functor for calculating the extended persistence hierarchy */
template <class Simplex> class ExtendedPersistenceHierarchy
{
public:
  using SimplicialComplex = topology::SimplicialComplex<Simplex>;
  using Vertex            = typename Simplex::VertexType;
  using SimplexPairing    = PersistencePairing<Vertex>;

  using EdgeType          = std::pair<Vertex, Vertex>;
  using Edges             = std::vector<EdgeType>;

  std::map<Simplex, Vertex> decomposeEdges( const SimplicialComplex& S )
  {
    std::set<Vertex> vertices;
    S.vertices( std::inserter( vertices, vertices.begin() ) );

    std::map<Vertex, Vertex>  criticalPointMapVertices;  // critical point map kept along with Union--Find structure
    std::map<Simplex, Vertex> criticalPointMapSimplices; // critical point map on simplex basis (return value)

    for( auto&& vertex : vertices )
      criticalPointMapVertices[vertex] = vertex;

    for( auto&& simplex : S )
    {
      if( simplex.dimension() != 1 )
        continue;

      Vertex u = *( simplex.begin() );
      Vertex v = *( simplex.begin() + 1 );

      auto youngerComponent = u;
      auto olderComponent   = v;

      if( youngerComponent == olderComponent )
      {
        criticalPointMapSimplices[ simplex ] = criticalPointMapVertices[ olderComponent ];
        continue;
      }
      else
      {
        auto index1 = S.index( Simplex(youngerComponent) );
        auto index2 = S.index( Simplex(olderComponent)   );

        if( index1 < index2 )
          std::swap( youngerComponent, olderComponent );

        Simplex creator = *( S.find( Simplex( youngerComponent ) ) );

        if( creator.data() == simplex.data() )
          criticalPointMapVertices[ youngerComponent ] = criticalPointMapVertices[ olderComponent ];

        criticalPointMapSimplices[ simplex ] = criticalPointMapVertices[ olderComponent ];
      }
    }

    return criticalPointMapSimplices;
  }

  std::pair<SimplicialComplex,
    topology::UnionFind<Vertex>
    //SDisjointSetForest<Vertex>
  >
    makeInterlevelSet( typename Simplex::DataType lower, typename Simplex::DataType upper,
                       const SimplicialComplex& S )
  {
    if( lower > upper )
      std::swap( lower, upper );

    std::vector<Simplex> simplices;

    std::copy_if( S.begin(), S.end(),
                  std::back_inserter( simplices ),
                  [&lower, &upper] ( const Simplex& s )
                  {
                    return s.data() >= lower && s.data() <= upper && s.dimension() <= 1;
                  } );

    SimplicialComplex K = SimplicialComplex( simplices.begin(),
                                             simplices.end() );

    std::set<Vertex> vertices;

    // FIXME: Need a better way of finding "proper" vertices in the simplicial
    // complex...
    for( auto&& simplex : simplices )
    {
      if( simplex.dimension() == 0 )
        vertices.insert( *simplex.begin() );
    }

    // FIXME: do I need this?
    //SDisjointSetForest<Vertex> uf( vertices.begin(), vertices.end() );
    topology::UnionFind<Vertex> uf( vertices.begin(), vertices.end()  );

    for( auto&& simplex : K )
    {
      if( simplex.dimension() == 1 )
      {
        Vertex u = *( simplex.begin() );
        Vertex v = *( simplex.begin() + 1 );

        if( !uf.contains( u ) || !uf.contains( v ) )
          continue;

        auto youngerComponent = uf.find( u );
        auto olderComponent   = uf.find( v );

        if( youngerComponent == olderComponent )
          continue;

        auto index1 = std::distance( S.begin(),
                                     S.find( Simplex( youngerComponent ) ) );

        auto index2 = std::distance( S.begin(),
                                     S.find( Simplex( olderComponent ) ) );

        if( index1 < index2 )
          std::swap( youngerComponent, olderComponent );

        uf.merge( youngerComponent, olderComponent );
      }
    }

    return std::make_pair( K, uf );
  }

  Edges operator()( const SimplicialComplex& simplicialComplex )
  {
    using namespace detail;

    // Extract {0,1}-simplices -----------------------------------------

    // Note that there is a range predicate for the simplicial complex class
    // that does essentially the same. However, the predicate is not stable
    // with respect to the filtration of the simplicial complex. Thus, to
    // extract the desired simplices, the internal function cannot be used.

    std::vector<Simplex> simplices;

    std::copy_if( simplicialComplex.begin(), simplicialComplex.end(),
                  std::back_inserter( simplices ),
                  [] ( const Simplex& s ) { return s.dimension() <= 1; } );

    SimplicialComplex S = SimplicialComplex( simplices.begin(),
                                             simplices.end() );

    // Persistence calculation -----------------------------------------

    // Get number of vertices stored in simplicial complex
    std::set<Vertex> vertices;
    S.vertices( std::inserter( vertices,
                               vertices.begin() ) );

    SimplexPairing pairing;

    Edges edges;

    auto criticalPointMap = decomposeEdges( S );

    // FIXME: is this map really required? The information given by the
    // critical point map seems to be the same...
    std::map<Vertex, Vertex> criticalPointsUF;
    for( auto&& vertex : vertices )
      criticalPointsUF[vertex] = vertex;

    // Required in order to obtain persistence pairs along with the
    // edges of the persistence hierarchy.
    topology::UnionFind<Vertex> uf( vertices.begin(), vertices.end() );

    for( auto&& simplex : S )
    {
      // Only edges can destroy a component
      if( simplex.dimension() != 1 )
        continue;

      auto u = *( simplex.begin() );
      auto v = *( simplex.begin() + 1 );

      // ---------------------------------------------------------------
      //
      // Ensure that the younger component is _always_ the first
      // component. A component is younger if its representative
      // vertex precedes the other vertex in the filtration.
      auto youngerComponent = uf.find( u );
      auto olderComponent   = uf.find( v );

      // If the component has already been merged by some other edge, we are
      // not interested in it any longer.
      if( youngerComponent == olderComponent )
        continue;

      {
        auto index1 = S.index( Simplex(youngerComponent) );
        auto index2 = S.index( Simplex(olderComponent)   );

        // The younger component has the _larger_ index as it is born _later_
        // in the filtration.
        if( index1 < index2 )
          std::swap( youngerComponent, olderComponent );
      }

      // Prepare information about creators ----------------------------
      //
      // Creator simplex for the simplex pairing below. I know that this
      // simplex must exist in the complex so I don't check for iterator
      // validity here.
      auto youngerCreator         = *( S.find( Simplex( youngerComponent ) ) );
      auto olderCreator           = *( S.find( Simplex( olderComponent ) ) );
      auto youngerCriticalSimplex = *( S.find( Simplex( criticalPointsUF[youngerComponent] ) ) );
      auto olderCriticalSimplex   = *( S.find( Simplex( criticalPointsUF[olderComponent]   ) ) );

      // Zero-persistence information; assign critical point of the
      // older component directly. This ensures that we are able to
      // obtain a proper decomposition.
      if( youngerCreator.data() == simplex.data() )
        criticalPointsUF[youngerComponent] = olderComponent;
      else
      {
        // Ensures that the oldest, highest/lowest critical simplex is
        // being used to calculate the interlevel set. Else, it may be
        // impossible for a critical point to be reached.
        if( S.index( youngerCriticalSimplex ) < S.index( olderCriticalSimplex ) )
          std::swap( youngerCriticalSimplex, olderCriticalSimplex );

        // At least one of the critical points of the two connected
        // components is trivial.
        if( !( criticalPointsUF[youngerComponent] != youngerComponent && criticalPointsUF[olderComponent] != olderComponent ) )
        {
          auto clsPair
            = makeInterlevelSet(
                olderCriticalSimplex.data(), simplex.data(),
                S );

          bool inSameComponent
            =   ( clsPair.second.contains( *olderCriticalSimplex.begin() ) && clsPair.second.contains( *youngerCriticalSimplex.begin() ) )
             && ( clsPair.second.find( *olderCriticalSimplex.begin() ) == clsPair.second.find( *youngerCriticalSimplex.begin() ) );

          if( inSameComponent )
          {
            using V = boost::graph_traits<AdjacencyGraph>::vertex_descriptor;

            boost::bimap<Vertex, SizeType> vim;
            AdjacencyGraph G;

            std::tie( vim, G )
              = extractZeroDimensionalAdjacencyGraph( clsPair.first );

            std::vector<V> p( boost::num_vertices( G ) );

            auto u = vim.left.at( *olderCriticalSimplex.begin() );
            auto v = vim.left.at( *youngerCriticalSimplex.begin() );
            p[u]   = u;

            boost::breadth_first_search( G,
                                         boost::vertex( u, G ),
                                         boost::visitor(
                                          boost::make_bfs_visitor(
                                            boost::record_predecessors( &p[0], boost::on_tree_edge() ) ) ) );

            std::set<Vertex> criticalPoints;

            while( true )
            {
              auto parent = p.at( v );
              auto s      = *S.find( Simplex( {vim.right.at(v),vim.right.at(parent)} ) );

              if( parent == v )
                break;

              v = parent;

              criticalPoints.insert( criticalPointMap.at( s ) );
            }

            // Exactly two critical points (i.e. the ones we were
            // looking for); hence, insert younger component as a
            // child of the youngest critical point.
            if( criticalPoints.size() == 2 )
            {
              edges.push_back( std::make_pair(
                criticalPointsUF[olderComponent],
                youngerComponent )
              );
            }

            // More critical points; connect the critical points
            // according to the usual persistence hierarchy.
            else
            {
              edges.push_back( std::make_pair(
                olderComponent,
                youngerComponent )
              );
            }
          }

          // Not in the same component; connect the critical points
          // according to the usual persistence hierarchy.
          else
          {
            edges.push_back( std::make_pair(
              olderComponent,
              youngerComponent )
            );
          }

        }

        // Case 3: Two critical points converge in one branch...
        else
        {
          throw std::runtime_error( "Not yet implemented" );

          // This code is old & dangerously untested. I am not even sure
          // whether this case can occur in data sets...
          #if 0
          std::cout << "- Critical simplex 1 (younger): " << youngerCriticalSimplex << "\n"
                    << "- Critical simplex 2 (older)  : " << olderCriticalSimplex << "\n";

          auto clsPair
            = makeInterlevelSet(
                olderCriticalSimplex.data(), simplex.data(),
                S );

          bool inSameComponent
            =  ( clsPair.second.contains( *olderCriticalSimplex.begin() ) && clsPair.second.contains( *youngerCriticalSimplex.begin() ) )
            && ( clsPair.second.find( *olderCriticalSimplex.begin() ) == clsPair.second.find( *youngerCriticalSimplex.begin() ) );

          //auto inSameComponent
          //  = clsPair.second.inSameComponent( *critical2.begin(), *critical1.begin() );

          std::cout << "Connected?\n"
                    << "   " << inSameComponent << "\n";

          if( !inSameComponent || inSameComponent )
          {
            edges.push_back( std::make_pair(
              olderComponent,
              youngerComponent )
            );

            std::cout << "# Edge: " << youngerComponent << " --> " << olderComponent << "\n";
          }
          #endif
        }

        // The youngest critical point along the current connected
        // component has been changed.
        criticalPointsUF[olderComponent] = youngerComponent;
      }

      // Store information in simplex pairing ------------------------

      #if 0
      pairing.add( Pair( creator,
                         cascade,
                         simplex ) );
      #endif

      // Actual merge ------------------------------------------------

      uf.merge( youngerComponent,
                olderComponent );
    }

    std::set<Vertex> roots;
    uf.roots( std::inserter( roots, roots.begin() ) );

    for( auto&& root : roots )
    {
      std::cout << "TRAVERSING ROOT " << root << "...\n"
                << "  CRITICAL POINT: " << criticalPointsUF[root] << "\n"
                << "  EDGE          : " << root << " <-- " << criticalPointsUF[root] << "\n";

      // TODO: What about this?
      //edges.push_back( std::make_pair(
      //  root,
      //  uf.getCritical( root ) )
      //);

      #if 0
      pairing.add( Pair( *S.find( root ) ) );
      #endif
    }

    //return std::make_pair( pairing, edges );
    return edges;
  }
};

} // namespace aleph

#endif
