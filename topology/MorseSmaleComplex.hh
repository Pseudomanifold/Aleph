#ifndef ALEPH_TOPOLOGY_MORSE_SMALE_COMPLEX__
#define ALEPH_TOPOLOGY_MORSE_SMALE_COMPLEX__

// TODO: Remove after debugging
#include <iostream>

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

      std::cerr << "No. higher neighbours: " << higherNeigbours.size() << "\n"
                << "No. lower  neighbours: " << lowerNeighbours.size() << "\n";

      std::cerr << "[" << vertex << "]: ";

      if( higherNeigbours.empty() )
        std::cerr << "Maximum\n";
      else if( lowerNeighbours.empty() )
        std::cerr << "Minimum\n";
      else
        std::cerr << "Saddle\n";
    }
  }
};

} // namespace topology

} // namespace aleph

#endif
