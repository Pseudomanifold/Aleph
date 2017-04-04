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
    auto n = M.vertices();

    for( decltype(n) i = 0; i < n; i++ )
    {
      auto v = M.vertex(i);

      auto higherNeigbours = M.getHigherNeighbours( *v );
      auto lowerNeighbours = M.getLowerNeighbours( *v );

      std::cerr << "No. higher neighbours: " << higherNeigbours.size() << "\n"
                << "No. lower  neighbours: " << lowerNeighbours.size() << "\n";

      std::cerr << "[" << i << "]: ";

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
