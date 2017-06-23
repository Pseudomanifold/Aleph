#ifndef ALEPH_TOPOLOGY_MESH_HH__
#define ALEPH_TOPOLOGY_MESH_HH__

#include <cassert>

#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <aleph/topology/UnionFind.hh>

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

    /**
      Collects all vertices of the given face. Vertex IDs will be
      returned in the order in which they are traversed along the
      face.
    */

    std::vector<Index> vertices() const
    {
      std::vector<Index> v;

      auto e = edge;

      do
      {
        v.push_back( e->source()->id );
        e = e->next;
      }
      while( e != edge );

      return v;
    }
  };

  struct Vertex
  {
    Index    id   = Index();
    Position x    = Position();
    Position y    = Position();
    Position z    = Position();
    Data     data = Data();

    HalfEdgePointer edge;
  };

  // Mesh attributes ---------------------------------------------------

  std::vector<Index> vertices() const noexcept
  {
    std::vector<Index> result;
    result.reserve( _vertices.size() );

    for( auto&& pair : _vertices )
      result.push_back( pair.first );

    return result;
  }

  std::size_t numVertices() const noexcept
  {
    return _vertices.size();
  }

  // TODO: Better return type?
  std::vector< std::vector<Index> > faces() const noexcept
  {
    std::unordered_set<FacePointer> faces;

    for( auto&& pair : _vertices )
    {
      auto&& vertex = pair.second;
      auto&& edge   = vertex->edge;
      auto&& face   = edge->face;

      if( face )
        faces.insert( face );
    }

    std::vector< std::vector<Index> > results;
    results.reserve( faces.size() );

    for( auto&& face : faces )
      results.push_back( face->vertices() );

    return results;
  }

  std::size_t numFaces() const noexcept
  {
    std::unordered_set<FacePointer> faces;

    for( auto&& pair : _vertices )
    {
      auto&& vertex = pair.second;
      auto&& edge   = vertex->edge;
      auto&& face   = edge->face;

      if( face )
        faces.insert( face );
    }

    return faces.size();
  }

  // Mesh modification -------------------------------------------------

  /** Adds a new vertex to the mesh */
  void addVertex( Position x, Position y, Position z, Data data = Data(), Index id = std::numeric_limits<Index>::max() )
  {
    Vertex v;

    v.id   = id == std::numeric_limits<Index>::max() ? std::max( _vertices.size(), _largestVertexID ) : id;
    v.x    = x;
    v.y    = y;
    v.z    = z;
    v.data = data;

    auto pair = _vertices.insert( std::make_pair( v.id,
                                  std::make_shared<Vertex>( v ) ) );

    if( !pair.second )
      throw std::runtime_error( "Vertex ID must be unique" );

    _largestVertexID = std::max( _largestVertexID, v.id );
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
    edges.reserve( static_cast<std::size_t>( std::distance( begin, end ) ) );

    for( InputIterator it = begin; it != end; ++it )
    {
      auto curr = it;
      auto next = std::next( it );

      if( next == end )
        next = begin;

      auto source = _vertices.at( *curr );   // Edge source vertex
      auto target = _vertices.at( *next );   // Edge target vertex
      auto edge   = getEdge( *curr, *next ); // Edge

      // Case 1: A new edge. Create a new edge and a new pair. Set edges
      // of source and target vertex correctly. Moreover, initialize the
      // paired edge with sensible default values.
      if( !edge )
      {
        edge      = std::make_shared<HalfEdge>();
        auto pair = std::make_shared<HalfEdge>();

        pair->vertex = source; // This is flipped by design: We point to the
                               // target vertex of the flipped edge. This is
                               // just the source vertex again.

        pair->pair   = edge;
        edge->pair   = pair;

        if( !source->edge )
          source->edge = edge;

        if( !target->edge )
          target->edge = pair;
      }

      // Ensures that the first edge that is specified for the new face
      // will be set as the outgoing edge of the face. This is not just
      // a 'cosmetic' choice but also ensures that vertex IDs for every
      // face are reported in the original order.
      if( !face->edge )
        face->edge = edge;

      assert( !edge->face );
      assert(  edge->pair );

      edge->face   = face;
      edge->vertex = target;

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

      edge->next = *next;
      edge->prev = *prev;
    }

    // Extend boundary -------------------------------------------------
    //
    // Traverse all vertices whose paired edges have an empty face. Any
    // of these edges is part of the boundary face.

    for( auto&& pair : _vertices )
    {
      auto&& vertex = pair.second;

      if( !vertex->edge || vertex->edge->pair->face )
        continue;

      auto curr = vertex->edge->target();
      auto edge = vertex->edge->pair;

      do
      {
        assert( !edge->face );
        auto edges = this->getEdges( *curr );

        for( auto&& e : edges )
        {
          if( !e->pair->face )
          {
            e->pair->next = edge;
            edge->prev    = e->pair;

            edge          = e->pair;
            curr          = e->pair->pair->target();
            break;
          }
        }
      }
      while( curr != vertex );

      // Close the loop around the boundary face by adding a pointer to
      // the identified edge.
      vertex->edge->pair->next = edge;
      edge->prev               = vertex->edge->pair;
    }
  }

  // Mesh queries ------------------------------------------------------

  /** Returns data stored at a certain vertex */
  Data data( Index id ) const
  {
    auto v = this->getVertex( id );
    return v->data;
  }

  /**
    The star of a vertex is defined as the mesh that contains all the
    triangles and edges of which the vertex is a face.
  */

  Mesh star( Index id ) const noexcept
  {
    Mesh M;

    auto vertex = this->getVertex( id );
    auto faces  = this->getFaces( *vertex );

    {
      std::unordered_set<Index> vertices;

      for( auto&& f : faces )
      {
        auto&& v = f->vertices();
        vertices.insert( v.begin(), v.end() );
      }

      for( auto&& v : vertices )
      {
        auto oldVertex = this->getVertex( v );

        M.addVertex( oldVertex->x, oldVertex->y, oldVertex->z,
                     oldVertex->data,
                     oldVertex->id );
      }
    }

    for( auto&& f : faces )
    {
      auto&& vertices = f->vertices();
      M.addFace( vertices.begin(), vertices.end() );
    }

    return M;
  }

  /**
    The link of a vertex is defined as all simplices in the closed star
    that are disjoint from the vertex. For 2-manifolds, this will yield
    a cycle of edges and vertices.

    This function will represent the cycle by returning all vertex IDs,
    in an order that is consistent with the orientation of the mesh.
  */

  std::vector<Index> link( Index id ) const noexcept
  {
    auto&& vertex     = this->getVertex( id );
    auto&& neighbours = this->getNeighbours( *vertex );

    std::vector<Index> result;
    result.reserve( neighbours.size() );

    for( auto&& neighbour : neighbours )
      result.push_back( neighbour->id );

    return result;
  }

  std::vector<Index> getLowerNeighbours( Index id ) const noexcept
  {
    auto&& vertex     = this->getVertex( id );
    auto&& neighbours = this->getNeighbours( *vertex );
    auto&& data       = vertex->data;

    neighbours.erase( std::remove_if( neighbours.begin(), neighbours.end(),
                                      [&data] ( const VertexPointer& neighbour )
                                      {
                                        return neighbour->data >= data;
                                      } ),
                      neighbours.end() );

    std::vector<Index> result;
    result.reserve( neighbours.size() );

    std::transform( neighbours.begin(), neighbours.end(), std::back_inserter( result ),
                    [] ( const VertexPointer& neighbour )
                    {
                      return neighbour->id;
                    } );

    return result;
  }

  std::vector<Index> getHigherNeighbours( Index id ) const noexcept
  {
    auto&& vertex     = this->getVertex( id );
    auto&& neighbours = this->getNeighbours( *vertex );
    auto&& data       = vertex->data;

    neighbours.erase( std::remove_if( neighbours.begin(), neighbours.end(),
                                      [&data] ( const VertexPointer& neighbour )
                                      {
                                        return neighbour->data <= data;
                                      } ),
                      neighbours.end() );

    std::vector<Index> result;
    result.reserve( neighbours.size() );

    std::transform( neighbours.begin(), neighbours.end(), std::back_inserter( result ),
                    [] ( const VertexPointer& neighbour )
                    {
                      return neighbour->id;
                    } );


    return result;
  }

  /**
    Checks whether an edge between two vertices that are identified by
    their index exists.
  */

  bool hasEdge( Index u, Index v ) const
  {
    auto&& source = this->getVertex(u);
    auto&& target = this->getVertex(v);

    auto neighbours = this->getNeighbours( *source );

    return std::find( neighbours.begin(), neighbours.end(), target ) != neighbours.end();
  }

  /** Counts the number of connected components */
  std::size_t numConnectedComponents() const noexcept
  {
    auto&& vertices = this->vertices();

    UnionFind<Index> uf( vertices.begin(), vertices.end() );

    for( auto&& vertex : vertices )
    {
      auto&& v      = _vertices.at( vertex );
      auto&& e      = v->edge;
      auto&& source = e->source()->id;
      auto&& target = e->target()->id;

      // Actually, the order does not matter one bit here...
      uf.merge( target, source );
    }

    std::vector<Index> roots;
    uf.roots( std::back_inserter( roots ) );

    return roots.size();
  }

private:

  /** Gets all vertices that are adjacent to a given vertex */
  std::vector<VertexPointer> getNeighbours( const Vertex& v ) const noexcept
  {
    std::vector<VertexPointer> neighbours;

    auto edge = v.edge;
    do
    {
      if( edge )
        neighbours.push_back( edge->target() );
      else
        break;

      edge = edge->prev->pair;
    }
    while( edge != v.edge && edge->source()->id == v.id );

    return neighbours;
  }

  /** Gets all edges that are incident on a given vertex. */
  std::vector<HalfEdgePointer> getEdges( const Vertex& v ) const noexcept
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

  /** Gets all faces that are incident on a given vertex */
  std::vector<FacePointer> getFaces( const Vertex& v ) const noexcept
  {
    std::vector<FacePointer> faces;

    auto edge = v.edge;
    do
    {
      if( edge && edge->face )
        faces.push_back( edge->face );
      else
        break;

      edge = edge->pair->next;
    }
    while( edge != v.edge );

    return faces;
  }

  /**
    Check whether a given (directed) edge already exists. If so,
    a pointer to the edge is being returned.
  */

  HalfEdgePointer getEdge( Index u, Index v ) const noexcept
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
    Returns vertex with the given ID. This is an internal function
    that is used to obtain all information about a vertex.
  */

  VertexPointer getVertex( Index id ) const
  {
    return _vertices.at( id );
  }

  /**
   Stores largest vertex ID. This is required in order to ensure
   that vertex IDs are not assigned multiple times when the user
   adds vertices one after the other.
  */

  Index _largestVertexID = Index();

  /**
    Stores all vertex pointers. This is sufficient to store the
    complete mesh.
  */

  std::unordered_map<Index, VertexPointer> _vertices;
};

} // namespace topology

} // namespace aleph

#endif
