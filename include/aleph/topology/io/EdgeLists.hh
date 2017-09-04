#ifndef ALEPH_TOPOLOGY_IO_EDGE_LISTS_HH__
#define ALEPH_TOPOLOGY_IO_EDGE_LISTS_HH__

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <aleph/utilities/String.hh>

namespace aleph
{

namespace topology
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

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K );
  }

  template <class SimplicialComplex> void operator()( std::ifstream& in, SimplicialComplex& K )
  {
    using namespace utilities;

    using Simplex           = typename SimplicialComplex::ValueType;
    using DataType          = typename Simplex::DataType;
    using VertexType        = typename Simplex::VertexType;

    std::string line;

    std::set<Simplex> vertices;
    std::vector<Simplex> edges;

    std::size_t lastID = 0;

    while( in )
    {
      std::getline( in, line );

      if( _trimLines )
        line = trim( line );

      // TODO: Make this configurable and permit splitting by different
      // tokens such as commas
      auto tokens = split( line );

      // Skip empty lines and comments
      if( line.empty() || std::find( _commentTokens.begin(), _commentTokens.end(), line.front() ) != _commentTokens.end() )
        continue;

      if( tokens.size() >= 2 )
      {
        VertexType u = VertexType();
        VertexType v = VertexType();

        // TODO: Make order of vertices & weights configurable?
        if( tokens[0].find_first_not_of( "0123456789" ) == std::string::npos )
        {
          u = convert<VertexType>( tokens[0] );
          v = convert<VertexType>( tokens[1] );
        }
        else
        {
          auto&& us = tokens[0];
          auto&& vs = tokens[1];

          if( _nodeLabels.find(us) == _nodeLabels.end() )
            _nodeLabels[us] = lastID++;

          if( _nodeLabels.find(vs) == _nodeLabels.end() )
            _nodeLabels[vs] = lastID++;

          u = static_cast<VertexType>( _nodeLabels.at(us) );
          v = static_cast<VertexType>( _nodeLabels.at(vs) );
        }

        DataType w = DataType();
        if( tokens.size() >= 3 && _readWeights )
          w = convert<DataType>( tokens[2] );

        edges.push_back( Simplex( { u, v }, w ) );

        vertices.insert( Simplex( u ) );
        vertices.insert( Simplex( v ) );
      }
      else
      {
        // TODO: Throw error?
      }
    }

    // Using a set has the advantage of ensuring that duplicate
    // simplices are deleted automatically. A duplicate simplex
    // is usually created by the input data set itself. It must
    // not be considered for any subsequent analysis.
    std::set<Simplex> simplices;

    simplices.insert( vertices.begin(), vertices.end() );
    simplices.insert( edges.begin(), edges.end() );

    K = SimplicialComplex( simplices.begin(), simplices.end() );
  }

private:
  std::vector<char> _commentTokens = { '#', '%', '\"', '*' };

  /**
    Optional set of node labels. These are read automatically in case an
    input data set does not use numeric labels. The idea is to map an ID
    in the form of a string to an index in corresponding graph.
  */

  std::map<std::string, std::size_t> _nodeLabels;

  bool _readWeights              = true;
  bool _trimLines                = true;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
