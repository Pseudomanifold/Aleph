#ifndef ALEPH_TOPOLOGY_MORSE_SMALE_COMPLEX__
#define ALEPH_TOPOLOGY_MORSE_SMALE_COMPLEX__

#include <algorithm>
#include <vector>

// TODO: Remove after debugging
#include <iostream>

#include <aleph/topology/UnionFind.hh>

namespace aleph
{

namespace topology
{

template <class Mesh> class MorseSmaleComplex
{
public:
  void operator()( const Mesh& M )
  {
    std::cerr << "\n";

    auto&& vertices = M.vertices();

    for( auto&& vertex : vertices )
    {
      auto higherNeigbours = M.getHigherNeighbours( vertex );
      auto lowerNeighbours = M.getLowerNeighbours( vertex );

      std::size_t nl = 0;
      std::size_t nu = 0;

      std::tie( nl, nu ) = contiguousSegments( M, vertex );

      std::cerr << "[" << vertex << "]: ";

      if( higherNeigbours.empty() )
        std::cerr << "Maximum\n";
      else if( lowerNeighbours.empty() )
        std::cerr << "Minimum\n";
      else
      {
        std::cerr << "Saddle\n";
        std::cerr << "  +: ";
        for( auto&& n : higherNeigbours )
          std::cerr << n << " " ;
        std::cerr << "\n";

        std::cerr << "  -: ";
        for( auto&& n : lowerNeighbours )
          std::cerr << n << " " ;
        std::cerr << "\n";
      }

      std::cerr << "  nl: " << nl << "\n"
                << "  nu: " << nu << "\n\n";

    }
  }
private:

  /**
    Calculates the number of contiguous segments in the link of
    a vertex, with respect to the given mesh.
  */

  static std::pair<std::size_t, std::size_t> contiguousSegments( const Mesh& M, typename Mesh::Index id )
  {
    using Index = typename Mesh::Index;

    auto data = M.data(id);
    auto link = M.link(id);

    std::vector<Index> upperLink;
    std::vector<Index> lowerLink;

    std::copy_if( link.begin(), link.end(), std::back_inserter( upperLink ), [&data, &M] ( Index& u ) { return M.data(u) >= data; } );
    std::copy_if( link.begin(), link.end(), std::back_inserter( lowerLink ), [&data, &M] ( Index& u ) { return M.data(u) <= data; } );

    auto&& numConnectedComponents = [&M, &id] ( const std::vector<Index>& link )
    {
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
    };

    return std::make_pair( numConnectedComponents( lowerLink ), numConnectedComponents( upperLink ) );
  }
};

} // namespace topology

} // namespace aleph

#endif
