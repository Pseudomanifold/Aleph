#ifndef ALEPH_TOPOLOGY_MORSE_SMALE_COMPLEX__
#define ALEPH_TOPOLOGY_MORSE_SMALE_COMPLEX__

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
    }
  }
};

} // namespace topology

} // namespace aleph

#endif
