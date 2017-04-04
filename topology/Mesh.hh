#ifndef ALEPH_TOPOLOGY_MESH_HH__
#define ALEPH_TOPOLOGY_MESH_HH__

#include <cassert>

#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>

namespace aleph
{

namespace topology
{

/**
  @class Mesh
  @brief Half-edge mesh data structure

  This data structure is capable of representing two-dimensional piecewise
  linear manifolds. In order to speed up standard queries, this class uses
  a standard half-edge data structure.
*/

template <class Position = float, class Data = float> class Mesh
{
public:
  struct Face;
  struct HalfEdge;
  struct Vertex;

  using Index           = std::size_t;

  using FacePointer     = std::shared_ptr<Face>;
  using HalfEdgePointer = std::shared_ptr<HalfEdge>;
  using VertexPointer   = std::shared_ptr<Vertex>;

  struct HalfEdge
  {
    FacePointer   face;
    VertexPointer vertex;

    HalfEdgePointer next; // Next half-edge (counter-clockwise)
    HalfEdgePointer prev; // Previous half-edge
    HalfEdgePointer pair; // Opposite half-edge

    VertexPointer source() const noexcept
    {
      return pair->vertex;
    }

    VertexPointer target() const noexcept
    {
      return vertex;
    }
  };

  struct Face
  {
    HalfEdgePointer edge;
  };

  struct Vertex
  {
    Position x = Position();
    Position y = Position();
    Position z = Position();
    Data     d = Data();

    HalfEdgePointer edge;
  };

  // Mesh modification -------------------------------------------------

  /** Adds a new vertex to the mesh */
  void addVertex( Position x, Position y, Position z, Data d = Data() )
  {
    Vertex v;

    v.x = x;
    v.y = y;
    v.z = z;
    v.d = d;

    _vertices.push_back( std::make_shared<Vertex>( v ) );
  }

  /**
    Adds a new face to the mesh. This function expects a range of vertex IDs
    that make up the face. The vertices of the face need to sorted correctly
    in order for the orientation to be consistent.
  */

  template <class InputIterator> void addFace( InputIterator begin, InputIterator end )
  {
    FacePointer face = std::make_shared<Face>();

    // Stores all half-edges created (or found) by this function in the
    // order in which they belong to the face.
    std::vector<HalfEdgePointer> edges;
    edges.reserve( std::distance( begin, end ) );

    for( InputIterator it = begin; it != end; ++it )
    {
      auto curr = it;
      auto next = std::next( it );

      if( next == end )
        next = begin;

      auto source = _vertices.at( *curr );   // Edge source vertex
      auto target = _vertices.at( *next );   // Edge target vertex
      auto edge   = getEdge( *curr, *next ); // Edge

      if( !edge )
      {
        edge      = std::make_shared<HalfEdge>();
        auto pair = std::make_shared<HalfEdge>();

        edge->face   = face;
        edge->pair   = pair;
        edge->vertex = target;
        pair->vertex = source;
        pair->pair   = edge;

        if( !source->edge )
          source->edge = edge;

        if( !target->edge )
          target->edge = pair;
      }
      else
      {
        assert( !edge->face );
        edge->face = face;
      }

      edges.push_back( edge );
    }

    // Set 'next' and 'prev' pointers correctly ------------------------
    //
    // We first traverse all edges that bound the current face. Here, it
    // should be possible to traverse the face directly, so we require a
    // proper pointer in both directions.

    for( auto itEdge = edges.begin(); itEdge != edges.end(); ++itEdge )
    {
      auto curr = itEdge;
      auto prev = std::prev( curr );
      auto next = std::next( curr );

      if( curr == edges.begin() )
        prev = std::prev( edges.end() );

      if( next == edges.end() )
        next = edges.begin();

      auto&& edge = *itEdge;

      assert( !edge->next );
      assert( !edge->prev );

      edge->next = *next;
      edge->prev = *prev;
    }
  }

  // Mesh queries ------------------------------------------------------

  std::vector<VertexPointer> getLowerNeighbours( const Vertex& v )
  {
    auto neighbours = this->getNeighbours( v );
    auto d          = v.d;

    neighbours.erase( std::remove_if( neighbours.begin(), neighbours.end(),
                                      [&d] ( const VertexPointer& neighbour )
                                      {
                                        return neighbour->d >= d;
                                      } ),
                      neighbours.end() );

    return neighbours;
  }

  std::vector<VertexPointer> getHigherNeighbours( const Vertex& v )
  {
    auto neighbours = this->getNeighbours( v );
    auto d          = v.d;

    neighbours.erase( std::remove_if( neighbours.begin(), neighbours.end(),
                                      [&d] ( const VertexPointer& neighbour )
                                      {
                                        return neighbour->d <= d;
                                      } ),
                      neighbours.end() );

    return neighbours;
  }

private:

  /** Gets all vertices that are adjacent to a given vertex */
  std::vector<VertexPointer> getNeighbours( const Vertex& v )
  {
    std::vector<VertexPointer> neighbours;

    auto edge = v.edge;
    do
    {
      if( edge )
        neighbours.push_back( edge->target );
      else
        break;

      edge = edge->pair->next;
    }
    while( edge != v.edge );

    return neighbours;
  }

  /** Gets all edges that are incident on a given vertex. */
  std::vector<HalfEdgePointer> getEdges( const Vertex& v )
  {
    std::vector<HalfEdgePointer> edges;

    auto edge = v.edge;
    do
    {
      if( edge )
        edges.push_back( edge );
      else
        break;

      edge = edge->pair->next;
    }
    while( edge != v.edge );

    return edges;
  }

  /**
    Check whether a given (directed) edge already exists. If so,
    a pointer to the edge is being returned.
  */

  HalfEdgePointer getEdge( Index u, Index v )
  {
    auto source = _vertices.at(u);           // Edge source vertex
    auto target = _vertices.at(v);           // Edge target vertex
    auto edges  = this->getEdges( *source ); // Incident edges

    auto itEdge = std::find_if( edges.begin(), edges.end(),
                                [&source, &target] ( const HalfEdgePointer& edge )
                                {
                                  return edge->source() == source && edge->target() == target;
                                } );

    if( itEdge != edges.end() )
      return *itEdge;
    else
      return nullptr;
  }

  /**
    Stores all vertex pointers. This is sufficient to store the
    complete mesh.
  */

  std::vector<VertexPointer> _vertices;
};

} // namespace topology

} // namespace aleph

#endif
