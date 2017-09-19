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

/**
  @class MatrixReader
  @brief Reads matrices in text format

  This reader class is meant to read matrices in text format. The matrix
  is supposed to consist of values that are separated by spaces. Each of
  the lines in the data stream is supposed to represent a *row*, whereas
  each individual value in a line is supposed to represent a *column*.

  The number of rows and columns must not vary over the file. An *empty*
  line is permitted, though. Likewise, lines starting with `#` will just
  be ignored. An example of a 3-by-3 matrix follows:

  \code
  0 1 2
  3 4 5
  6 7 8
  \endcode

  This file format is particularly useful for representing grey-scale
  images in an easily-understandable format. The reader stores all of
  the information of the data and permits queries about the width and
  height of the resulting matrix. The expansion behaviour may also be
  modified: by default, two-dimensional simplices (triangles) will be
  added to the simplicial complex. This can be triggered by setting a
  flag via `MatrixReader::addTriangles()`.
*/

class MatrixReader
{
public:

  /**
    Reads a simplicial complex from a file, using the default maximum
    functor for weight assignment. If you want to change the functor,
    please refer to the overloaded variant of this method.

    @param filename Input filename
    @param  K       Simplicial complex

    @see operator()( const std::string&, SimplicialComplex& )
  */

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;

    this->operator()( filename, K, [] ( DataType a, DataType b ) { return std::max(a,b); } );
  }

  /**
    Reads a simplicial complex from a file while supporting arbitrary
    functors for weight assignment. The functor needs to support this
    interface:

    \code{.cpp}
    using SimplexType = typename SimplicialComplex::ValueType;
    using DataType    = typename Simplex::DataType;

    DataType Functor::operator()( DataType a, DataType b )
    {
      // Do something with a and b. Typically, this would be calculating
      // either the minimum or the maximum...
      return std::max(a, b);
    }
    \endcode

    Please refer to the documentation of SimplicialComplexReader::operator()( const std::string& SimplicialComplex&, Functor )
    for more details.
  */

  template <class SimplicialComplex, class Functor> void operator()( const std::string& filename, SimplicialComplex& K, Functor f )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K, f );
  }

  /** @overload operator()( const std::string&, SimplicialComplex& ) */
  template <class SimplicialComplex> void operator()( std::istream& in, SimplicialComplex& K )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;

    this->operator()( in, K, [] ( DataType a, DataType b ) { return std::max(a,b); } );
  }

  /** @overload operator()( const std::string&, SimplicialComplex&, SimplicialComplex&, Functor ) */
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

    if( _addTriangles )
    {
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
    }

    K = SimplicialComplex( simplices.begin(), simplices.end() );
  }

  /** @returns Height of matrix that was read last */
  std::size_t height() const noexcept { return _height; }

  /** @returns Width of matrix that was read last */
  std::size_t width()  const noexcept { return _width;  }

  /**
    Configures expansion behaviour of the reader and determines whether
    triangles should be added or not.
  */

  void addTriangles( bool value = true )
  {
    _addTriangles = value;
  }

private:
  std::size_t _height = 0;
  std::size_t _width  = 0;

  bool _addTriangles  = true;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
