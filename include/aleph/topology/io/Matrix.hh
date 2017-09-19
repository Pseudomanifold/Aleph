#ifndef ALEPH_TOPOLOGY_IO_MATRIX_HH__
#define ALEPH_TOPOLOGY_IO_MATRIX_HH__

#include <algorithm>
#include <istream>
#include <fstream>
#include <string>
#include <vector>

#include <aleph/utilities/String.hh>

namespace aleph
{

namespace topology
{

namespace io
{

class MatrixReader
{
public:
  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;

    this->operator()( filename, K, [] ( DataType a, DataType b ) { return std::max(a,b); } );
  }

  template <class SimplicialComplex, class Functor> void operator()( const std::string& filename, SimplicialComplex& K, Functor f )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K, f );
  }

  template <class SimplicialComplex> void operator()( std::istream& in, SimplicialComplex& K )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;

    this->operator()( in, K, [] ( DataType a, DataType b ) { return std::max(a,b); } );
  }

  template <class SimplicialComplex, class Functor> void operator()( std::istream& in, SimplicialComplex& K, Functor f )
  {
    auto position = in.tellg();

    std::size_t height = 0;
    std::size_t width  = 0;

    using namespace aleph::utilities;

    {
      std::string line;
      while( std::getline( in, line ) )
      {
        line        = trim( line );
        auto tokens = split( line );

        if( width == 0 )
          width = tokens.size();
        else if( width != tokens.size() )
          throw std::runtime_error( "Format error: number of columns must not vary" );

        ++height;
      }

      in.clear();
      in.seekg( position );
    }

    _height = height;
    _width  = width;

    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;
    using VertexType = typename Simplex::VertexType;

    std::vector<DataType> values;
    values.reserve( _height * _width );

    std::copy( std::istream_iterator<DataType>( in ), std::istream_iterator<DataType>(),
               std::back_inserter( values ) );

    std::vector<Simplex> simplices;

    // Vertices --------------------------------------------------------

    {
      VertexType v = VertexType();

      for( auto&& value : values )
        simplices.push_back( Simplex( v++, value ) );
    }

    // Horizontal edges ------------------------------------------------

    for( std::size_t y = 0; y < _height; y++ )
    {
      for( std::size_t x = 0; x < _width - 1; x++ )
      {
        // current:   (x,  y) --> width * y + x
        // neighbour: (x+1,y) --> width * y + x+1
        auto u = static_cast<VertexType>( width * y + x   );
        auto v = static_cast<VertexType>( width * y + x+1 );
        auto w = f( values.at(u), values.at(v) );

        simplices.push_back( Simplex( {u,v}, w ) );
      }
    }

    // Vertical edges --------------------------------------------------

    for( std::size_t y = 0; y < _height - 1; y++ )
    {
      for( std::size_t x = 0; x < _width; x++ )
      {
        // current:   (x,y  ) --> width * (y  ) + x
        // neighbour: (x,y+1) --> width * (y+1) + x
        auto u = static_cast<VertexType>( width * (y  ) + x );
        auto v = static_cast<VertexType>( width * (y+1) + x );
        auto w = f( values.at(u), values.at(v) );

        simplices.push_back( Simplex( {u,v}, w ) );
      }
    }

    // Diagonal edges --------------------------------------------------

    for( std::size_t y = 0; y < _height - 1; y++ )
    {
      for( std::size_t x = 0; x < _width - 1; x++ )
      {
        // current:   (x  ,y  ) --> width * (y  ) + x
        // neighbour: (x+1,y+1) --> width * (y+1) + x+1
        auto u = static_cast<VertexType>( width * (y  ) + x   );
        auto v = static_cast<VertexType>( width * (y+1) + x+1 );
        auto w = f( values.at(u), values.at(v) );

        simplices.push_back( Simplex( {u,v}, w ) );
      }
    }

    // Triangles -------------------------------------------------------

    for( std::size_t y = 0; y < _height - 1; y++ )
    {
      for( std::size_t x = 0; x < _width - 1; x++ )
      {
        /*
          [a] (x,y  ) o---o (x+1,y  ) [b]
                      |\  |
                      | \ |
                      |  \|
          [d] (x,y+1) o---o (x+1,y+1) [c]

          The two triangles are formed by the respective corner points.

          [a] (x  ,y  ):  width * (y  ) + x
          [b] (x+1,y  ):  width * (y  ) + x+1
          [c] (x+1,y+1):  width * (y+1) + x+1
          [d] (x  ,y+1):  width * (y+1) + x
        */

        auto a = static_cast<VertexType>( width * (y  ) + x   );
        auto b = static_cast<VertexType>( width * (y  ) + x+1 );
        auto c = static_cast<VertexType>( width * (y+1) + x+1 );
        auto d = static_cast<VertexType>( width * (y+1) + x   );

        auto v = f( f( values.at(a), values.at(b) ), values.at(c) );
        auto w = f( f( values.at(a), values.at(c) ), values.at(d) );

        simplices.push_back( Simplex( {a,b,c}, v ) );
        simplices.push_back( Simplex( {a,c,d}, w ) );
      }
    }

    K = SimplicialComplex( simplices.begin(), simplices.end() );
  }

  std::size_t height() const noexcept { return _height; }
  std::size_t width()  const noexcept { return _width;  }

private:
  std::size_t _height = 0;
  std::size_t _width  = 0;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
