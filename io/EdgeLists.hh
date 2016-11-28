#ifndef ALEPH_IO_EDGE_LISTS_HH__
#define ALEPH_IO_EDGE_LISTS_HH__

#include <fstream>
#include <list>

#include "SimplicialComplex.hh"

namespace aleph
{

namespace io
{

class EdgeListReader
{
public:

  bool readWeights() const noexcept { return _readWeights; }
  bool trimLines()   const noexcept { return _trimLines;   }

  void setReadWeights( bool value = true ) noexcept { _readWeights = value; }
  void setTrimLines( bool value = true )   noexcept { _trimLines = value; }

  template <class DataType,
    class VertexType> SimplicialComplex< Simplex<DataType, VertexType> > operator()( std::ifstream& in )
  {
  }

private:
  std::list<char> _commentTokens = '#%';

  bool _readWeights              = true;
  bool _trimLines                = true;
};

}

#endif
