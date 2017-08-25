#ifndef ALEPH_TOPOLOGY_IO_HDF5_HH__
#define ALEPH_TOPOLOGY_IO_HDF5_HH__

// FIXME: guard this inclusion
#include <H5Cpp.h>

#include <algorithm>
#include <string>
#include <vector>

#include <cassert>

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class HDF5SimpleDataSpaceReader
  @brief Supports reading simple data spaces from HDF5 files

  This class provides a parser for extracting *simple data spaces* from
  HDF5 files. In other words, this class can extract scalar fields from
  HDF5 files. The class can only extract one field at a time. Moreover,
  it requires knowledge about which group and which data set to parse.
*/

class HDF5SimpleDataSpaceReader
{
public:

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    using Simplex  = typename SimplicialComplex::ValueType;
    using DataType = typename Simplex::DataType;

    this->operator()( filename, K, [] ( DataType a, DataType b ) { return std::max(a,b); } );
  }

  template <class SimplicialComplex, class Functor> void operator()( const std::string& filename, SimplicialComplex& K, Functor f )
  {
    using namespace H5;

    H5File file( filename, H5F_ACC_RDONLY );

    auto&& group       = file.openGroup( _groupName );
    auto&& dataSet     = group.openDataSet( _dataSetName );
    auto&& dataSpace   = dataSet.getSpace();
    auto&& dimension   = dataSpace.getSimpleExtentNdims();
    auto&& numVertices = dataSpace.getSimpleExtentNpoints();

    if( dimension != 2 )
      return;

    std::size_t width  = 0;
    std::size_t height = 0;

    {
      hsize_t* dimensions = new hsize_t[ dimension ];

      dataSpace.getSimpleExtentDims( dimensions,
                                     nullptr );

      width  = dimensions[0];
      height = dimensions[1];

      assert( width * height == std::size_t( numVertices ) );
    }

    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;
    using VertexType = typename Simplex::VertexType;

    auto n    = width * height;
    auto data = read<DataType>( dataSet, n );

    std::vector<Simplex> simplices;

    // 0-skeleton ------------------------------------------------------
    //
    // This is easy, as we can just copy all data directly into the
    // simplex data structure.

    for( std::size_t i = 0; i < n; i++ )
      simplices.push_back( Simplex( static_cast<VertexType>( i ), data.at(i) ) );

    // 1-skeleton & 2-skeleton -----------------------------------------
    //
    // This is a little bit more involved, but an auxiliary functor
    // helps us map coordinates to indices directly.

    auto coordinatesToIndex = [&width] ( std::size_t x, std::size_t y )
    {
      return y * width + x;
    };

    for( std::size_t y = 0; y < height; y++ )
    {
      for( std::size_t x = 0; x < width; x++ )
      {
        auto&& i = coordinatesToIndex( x, y );

        std::vector<std::size_t> neighbourCoordinates;
        std::vector<std::size_t> triangleCoordinates;

        if( x > 0 )
          neighbourCoordinates.push_back( coordinatesToIndex( x-1, y ) );
        if( x < width - 1 )
          neighbourCoordinates.push_back( coordinatesToIndex( x+1, y ) );
        if( y > 0 )
          neighbourCoordinates.push_back( coordinatesToIndex( x  , y-1 ) );
        if( y < height - 1 )
          neighbourCoordinates.push_back( coordinatesToIndex( x  , y+1 ) );

        if( x < width - 1 && y > 0 )
        {
          neighbourCoordinates.push_back( coordinatesToIndex( x+1, y-1 ) );

          triangleCoordinates.push_back( coordinatesToIndex( x,   y   ) );
          triangleCoordinates.push_back( coordinatesToIndex( x,   y-1 ) );
          triangleCoordinates.push_back( coordinatesToIndex( x+1, y   ) );
          triangleCoordinates.push_back( coordinatesToIndex( x+1, y-1 ) );
        }

        // 1-skeleton --------------------------------------------------
        //
        // The check below ensures that an edge is only added once. This
        // is faster than using e.g. std::set.

        for( auto&& n : neighbourCoordinates )
        {
          auto u = i;
          auto v = n;

          if( u > v )
            simplices.push_back( Simplex( {u,v}, f( data[u], data[v] ) ) );
        }

        if( !triangleCoordinates.empty() )
        {
          auto&& u = triangleCoordinates.at(0);
          auto&& v = triangleCoordinates.at(1);
          auto&& w = triangleCoordinates.at(2);
          auto&& x = triangleCoordinates.at(3);

          simplices.push_back( Simplex( {u,x,w}, f( data[u], f( data[x], data[w] ) ) ) );
          simplices.push_back( Simplex( {u,v,x}, f( data[u], f( data[v], data[x] ) ) ) );
        }
      }
    }

    K = SimplicialComplex( simplices.begin(), simplices.end() );
  }

  // Getters -----------------------------------------------------------

  std::string groupName() const noexcept   { return _groupName; }
  std::string dataSetName() const noexcept { return _dataSetName; }

  // Setters -----------------------------------------------------------

  void setGroupName( const std::string& name ) noexcept   { _groupName = name;   }
  void setDataSetName( const std::string& name ) noexcept { _dataSetName = name; }

private:

  /**
    Auxiliary function for reading data from an HDF5 file and *directly*
    into a vector. This function does not perform additional type checks
    but merely attempts a direct conversion.
  */

  template <class T> std::vector<T> read( const H5::DataSet& dataSet, std::size_t size )
  {
    std::vector<T> data;

    auto dataType = dataSet.getDataType();

    std::vector<H5::PredType> doubleTypes = { H5::PredType::NATIVE_DOUBLE,
                                              H5::PredType::IEEE_F64LE,
                                              H5::PredType::IEEE_F64BE };

    std::vector<H5::PredType> floatTypes  = { H5::PredType::NATIVE_FLOAT,
                                              H5::PredType::IEEE_F32LE,
                                              H5::PredType::IEEE_F32BE };

    if( std::find( doubleTypes.begin(), doubleTypes.end(), dataType ) != doubleTypes.end() )
    {
      std::vector<double> data_( size );

      dataSet.read( &data_[0], dataType );
      data.assign( data_.begin(), data_.end() );
    }
    else if( std::find( floatTypes.begin(), floatTypes.end(), dataType ) != floatTypes.end() )
    {
      std::vector<float> data_( size );

      dataSet.read( &data_[0], dataType );
      data.assign( data_.begin(), data_.end() );
    }
    else
    {
      // FIXME: handle unknown data type
    }

    return data;
  }

  std::string _groupName    = "/";
  std::string _dataSetName  = "YField";
};

} // namespace io

} // namespace topology

} // namespace aleph


#endif
