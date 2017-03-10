#ifndef ALEPH_TOPOLOGY_IO_PAJEK_HH__
#define ALEPH_TOPOLOGY_IO_PAJEK_HH__

#include <fstream>
#include <string>
#include <vector>

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class PajekReader
  @brief Parses files in Pajek format

  This is a simple reader for graphs in Pajek format. It supports loading
  Pajek files with vertex labels and edge weights.
*/

class PajekReader
{
public:

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K );
  }

  template <class SimplicialComplex> void operator()( std::ifstream& in, SimplicialComplex& K )
  {
    Mode mode = Mode::Unspecified;

    std::string line;
    while( std::getline( in, line ) )
    {
    }
  }

private:

  enum class Mode
  {
    Vertices,
    Edges,
    Unspecified
  };

  // Local storage -----------------------------------------------------
  //
  // Variables in this section contain the result of the last parsing
  // process. This is useful when clients are querying attributes.

  std::vector<std::string> _labels;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
