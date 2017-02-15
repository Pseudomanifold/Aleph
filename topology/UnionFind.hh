#ifndef ALEPH_TOPOLOGY_UNION_FIND_HH__
#define ALEPH_TOPOLOGY_UNION_FIND_HH__

#include <algorithm>

#include <unordered_map>
#include <unordered_set>

namespace aleph
{

namespace topology
{

template <class Vertex> class UnionFind
{
public:

  /**
    Creates a new Union--Find data structure and initializes a set of items
    with themselves.
  */

  template <class InputIterator> UnionFind( InputIterator begin,
                                            InputIterator end )
  {
    for( auto it = begin; it != end; ++it )
      _parent[ *it ]   = *it;
  }

  /**
    Merges a given vertex $u$ into the set corresponding to vertex $v$. Note
    that the merge is directional.
  */

  void merge( Vertex u, Vertex v ) noexcept
  {
    if( u != v )
      _parent[ this->find( u ) ] = this->find( v );
  }

  /**
    Finds the parent of a given vertex. Performs path collapse operations, if
    necessary. The function will throw if it encounters an unknown vertex.
  */

  Vertex find( Vertex u )
  {
    if( _parent.at(u) == u )
      return u;
    else
    {
      // Path collapse!
      _parent[u] = this->find( _parent.at( u ) );
      return _parent.at( u );
    }
  }

  /**
    Enumerates all roots, i.e. all sets that have themselves as a parent
    vertex, and stores it using an output iterator. Vertices will appear
    in random order.
  */

  template <class OutputIterator> void roots( OutputIterator result )
  {
    std::unordered_set<Vertex> r;

    for( auto&& pair : _parent )
    {
      auto&& u = pair.first;

      if( this->find( u ) == u )
        r.insert( u );
    }

    std::copy( r.begin(), r.end(), result );
  }

  /**
    Gets all vertices that make up a given connected component. The
    parameter of this function should be a root vertex, i.e. one of
    the creator vertices in the data structure.

    Vertices will appear in random order.
  */

  template <class OutputIterator> void get( Vertex v, OutputIterator result )
  {
    std::unordered_set<Vertex> r;

    auto&& parent = this->find(v);
    for( auto&& pair : _parent )
    {
      if( this->find( pair.first ) == parent )
        r.insert( pair.first );
    }

    std::copy( r.begin(), r.end(), result );
  }

private:

  /** Stores the usual parent--child relationship */
  std::unordered_map<Vertex, Vertex> _parent;
};

} // namespace topology

} // namespace aleph

#endif
