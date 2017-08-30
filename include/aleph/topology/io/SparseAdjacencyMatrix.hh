#ifndef ALEPH_TOPOLOGY_IO_SPARSE_ADJACENCY_MATRIX_HH__
#define ALEPH_TOPOLOGY_IO_SPARSE_ADJACENCY_MATRIX_HH__

#include <aleph/utilities/String.hh>

// FIXME: remove after debugging
#include <iostream>

#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>

namespace aleph
{

namespace topology
{

namespace io
{

class SparseAdjacencyMatrixReader
{
public:
  template <class SimplicialComplex> void operator()( const std::string& filenameMatrix,
                                                      const std::string& filenameIndicator,
                                                      std::vector<SimplicialComplex>& result )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using VertexType = typename Simplex::VertexType;
    using Edge       = std::pair<VertexType, VertexType>;

    std::vector<Edge> edges;

    {
      std::ifstream in( filenameMatrix );
      if( !in )
        throw std::runtime_error( "Unable to read input adjacency matrix file" );

      std::string line;
      while( std::getline( in, line ) )
      {
        using namespace aleph::utilities;

        auto tokens = split( line, _separator );

        if( tokens.size() == 2 )
        {
          auto u = convert<VertexType>( tokens.front() );
          auto v = convert<VertexType>( tokens.back()  );

          edges.push_back( std::make_pair(u, v) );
        }
        else
        {
          // TODO: throw error?
        }
      }

      std::cerr << __PRETTY_FUNCTION__ << ": Read " << edges.size() << " edges\n";
    }

    std::vector<std::string> graphs;

    {
      std::ifstream in( filenameIndicator );
      if( !in )
        throw std::runtime_error( "Unable to read input graph indicator file" );

      std::string line;
      while( std::getline( in, line ) )
      {
        using namespace aleph::utilities;
        line = trim( line );

        graphs.push_back( line );
      }

      std::cerr << __PRETTY_FUNCTION__ << ": Read " << graphs.size() << " graph indicators\n";
    }

    {
      std::set<std::string> G( graphs.begin(), graphs.end() );
      std::cerr << __PRETTY_FUNCTION__ << ": Read " << G.size() << " different graphs\n";
    }

    (void) result;
  }

private:

  // TODO: make configurable
  std::string _separator = ",";
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
