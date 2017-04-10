#ifndef ALEPH_TOPOLOGY_MORSE_SMALE_COMPLEX__
#define ALEPH_TOPOLOGY_MORSE_SMALE_COMPLEX__

// TODO: Remove after debugging
#include <iostream>

#include "topology/UnionFind.hh"

namespace aleph
{

namespace topology
{

template <class Mesh> class MorseSmaleComplex
{
public:
  void operator()( const Mesh& M )
  {
    auto&& vertices = M.vertices();

    for( auto&& vertex : vertices )
    {
      auto higherNeigbours = M.getHigherNeighbours( vertex );
      auto lowerNeighbours = M.getLowerNeighbours( vertex );

      std::cerr << "[" << vertex << "]: ";

      if( higherNeigbours.empty() )
      {
        std::cerr << "Maximum\n"
                  << "  k = " << contiguousSegments( M, vertex ) << "\n";
      }
      else if( lowerNeighbours.empty() )
        std::cerr << "Minimum\n"
                  << "  k = " << contiguousSegments( M, vertex ) << "\n";
      else
      {
        std::cerr << "Saddle\n"
                  << "  k = " << contiguousSegments( M, vertex ) << "\n";

        std::cerr << "  + ";
        for( auto&& neighbour : higherNeigbours )
          std::cerr << neighbour << " ";
        std::cerr << "\n";

        std::cerr << "  - ";
        for( auto&& neighbour : lowerNeighbours )
          std::cerr << neighbour << " ";
        std::cerr << "\n";
      }
    }
  }
private:

  /**
    Calculates the number of contiguous segments in the link of
    a vertex, with respect to the given mesh.
  */

  static std::size_t contiguousSegments( const Mesh& M, typename Mesh::Index id )
  {
    using Index = typename Mesh::Index;

    auto link = M.link(id);
    auto curr = link.begin();
    auto next = std::next( curr );

    UnionFind<Index> uf( link.begin(), link.end() );

    for( ; curr != link.end(); )
    {
      if( next == link.end() )
        next = link.begin();

      auto u    = id;
      auto v    = *curr++;
      auto w    = *next++;

      if( M.hasEdge(u,v) )
      {
        if( M.hasEdge(v,w) )
          uf.merge(v,w);
      }
      else
        throw std::runtime_error( "Edge to link centre vertex does not exist" );
    }

    std::vector<Index> roots;
    uf.roots( std::back_inserter(roots) );

    return roots.size();
  }
};

} // namespace topology

} // namespace aleph

#endif
