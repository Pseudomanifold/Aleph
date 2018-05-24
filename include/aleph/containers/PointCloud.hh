#ifndef ALEPH_CONTAINERS_POINT_CLOUD_HH__
#define ALEPH_CONTAINERS_POINT_CLOUD_HH__

#include <algorithm>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <string>

#include <cstddef>

#include <aleph/utilities/String.hh>

namespace aleph
{

namespace containers
{

template <class T> class PointCloud
{
public:

  // Exporting the type of elements stored in the point cloud. This is
  // used by other algorithms to prevent casting.
  using ElementType = T;

  // Ditto for the index type. It might make sense to leave this
  // configurable, but I do not see a pressing reason to do so.
  using IndexType = std::size_t;

  PointCloud()
    : _n( 0 )
    , _d( 0 )
    , _points( nullptr )
  {
  }

  PointCloud( std::size_t n, std::size_t d )
    : _n( n )
    , _d( d )
    , _points( new T[ _n * _d ] )
  {
    // Zero-initialization. This may not be the most efficient way of
    // doing it, in particular if a client has some data to pass, but
    // it ensures consistency.
    std::fill( _points, _points + _n * _d, T() );
  }

  PointCloud( const PointCloud& other )
    : _n( other._n )
    , _d( other._d )
    , _points( new T[ _n * _d ] )
  {
    std::copy( other._points, other._points + _n * _d, _points );
  }

  PointCloud( PointCloud&& other )
    : PointCloud()
  {
    swap( *this, other );
  }

  PointCloud& operator=( PointCloud other )
  {
    swap( *this, other );
    return *this;
  }

  ~PointCloud()
  {
    delete[] _points;
  }

  friend void swap( PointCloud& pc1, PointCloud& pc2 ) noexcept
  {
    using std::swap;

    swap( pc1._points, pc2._points );
    swap( pc1._n,      pc2._n );
    swap( pc1._d,      pc2._d );
  }

  // Equality comparison -----------------------------------------------

  bool operator==( const PointCloud<T>& other ) const noexcept
  {
    return    _n == other._n
           && _d == other._d
           && std::equal( _points, _points + _n * _d, other._points );
  }

  // Attributes --------------------------------------------------------

  IndexType size() const noexcept
  {
    return _n;
  }

  IndexType dimension() const noexcept
  {
    return _d;
  }

  bool empty() const noexcept
  {
    return _n == 0;
  }

  // Point access ------------------------------------------------------

  // This is slightly evil. The function is not really "bit-wise"
  // constant because it still permits modifying the pointer that
  // is being returned.
  T* data() const noexcept
  {
    return _points;
  }

  /**
    Sets $i$th point of point cloud. Throws if the number of dimensions
    does not match the number of dimensions in the point cloud.
  */

  template <class InputIterator> void set( IndexType i,
                                           InputIterator begin, InputIterator end )
  {
    if( i >= this->size() )
      throw std::runtime_error( "Invalid index" );

    auto distance = std::distance( begin, end );

    if( static_cast<IndexType>( distance ) != this->dimension() )
      throw std::runtime_error( "Incorrect number of dimensions" );

    std::copy( begin, end, _points + this->dimension() * i );
  }

  /** @overload set() */
  void set( IndexType i, const std::initializer_list<T>& il )
  {
    this->set( i, il.begin(), il.end() );
  }

  /**
    Gets the $i$th point of the point cloud. It will be stored via an
    output iterator. Incorrect indices will result in an exception.
  */

  template <class OutputIterator> void get( IndexType i,
                                            OutputIterator result ) const
  {
    if( i >= this->size() )
      throw std::runtime_error( "Invalid index" );

    auto offset = i * this->dimension();

    std::copy( _points + offset, _points + offset + this->dimension(),
               result );
  }

  /**
    Stores the $i$th point of the point cloud in a vector and returns
    said vector. This interface is not very efficient but is provided
    as a fallback/simplified variant.
  */

  std::vector<T> operator[]( IndexType i ) const
  {
    std::vector<T> p;

    this->get( i,
               std::back_inserter( p ) );

    return p;
  }

  // Operations --------------------------------------------------------

  /**
    Permits the concatenation of two point clouds. To this end, all of
    the points in the first point cloud will be appended to the points
    in the second point cloud, thereby forming a new point cloud. This
    new point cloud is then returned. Note that the dimensions of both
    point clouds have to be equal. Otherwise, an error is thrown.
  */

  PointCloud operator+( const PointCloud& other ) const
  {
    if( this->dimension() != other.dimension() )
      throw std::runtime_error( "The dimensions of both point clouds have to coincide" );

    auto d = this->dimension();
    auto n = this->size() + other.size();

    PointCloud result(n, d);

    decltype(n) i = 0; // total index with respect to result point cloud
    decltype(n) j = 0; // local index with respect to current point cloud

    for( j = 0; j < this->size(); j++, i++ )
    {
      auto&& p = this->operator[](j);
      result.set( i, p.begin(), p.end() );
    }

    for( j = 0; j < other.size(); j++, i++ )
    {
      auto&& p = other[j];
      result.set(i, p.begin(), p.end() );
    }

    return result;
  }

private:
  std::size_t _n; ///< Number of points
  std::size_t _d; ///< Dimension

  T* _points;
};

/**
  Loads a new point cloud from a file. The file is supposed to be in
  ASCII format. Each row must specify one item of the data set.  The
  different attributes of each item are assumed to be separated by a
  comma or white-space characters.
*/

template<class T> PointCloud<T> load( const std::string& filename )
{
  std::ifstream in( filename );

  if( !in )
    return PointCloud<T>();

  auto lines = std::count( std::istreambuf_iterator<char>( in ),
                           std::istreambuf_iterator<char>(),
                           '\n' );

  in.clear();
  in.seekg( 0 );

  if( lines <= 0 )
    return PointCloud<T>();

  // Recount all lines in the file, removing comments and empty lines
  // from the calculation. Else, the code below allocates many points
  // that will remain empty.
  {
    std::string line;
    while( std::getline( in, line ) )
    {
      line = utilities::trim( line );
      if( line.empty() || line.front() == '#' )
        lines -= 1;
    }
  }

  in.clear();
  in.seekg( 0 );

  std::size_t i = 0;
  std::size_t d = 0;
  std::size_t n = static_cast<std::size_t>( lines );

  std::string line;

  PointCloud<T> pointCloud;

  while( std::getline( in, line ) )
  {
    line        = utilities::trim( line );
    auto tokens = utilities::split( line, std::string( "[:;,[:space:]]+" ) );

    // Skip comment lines or empty lines; while this is somewhat
    // superfluous in most files (at least it is unlikely that a
    // line in the middle of the file will be empty), the loader
    // should handle empty lines at the end of the file.
    if( line.empty() || line.front() == '#' || tokens.empty() )
      continue;

    if( d == 0 )
    {
      d          = tokens.size();
      pointCloud = PointCloud<T>( n, d );
    }

    std::vector<T> coordinates;
    coordinates.reserve( d );

    for( auto&& token : tokens )
    {
      T coordinate = utilities::convert<T>( token );
      coordinates.push_back( coordinate );
    }

    pointCloud.set( i,
                    coordinates.begin(), coordinates.end() );

    ++i;
  }

  return pointCloud;
}

} // namespace containers

} // namespace aleph

/**
  Output operator for writing a point cloud to an `std::ostream`. The
  Individual attributes of the point cloud are separated by tabs, and
  each point is on a single line.

  @param o          Output stream
  @param pointCloud Point cloud to store

  @returns Modified output stream
*/

template <class T> std::ostream& operator<<( std::ostream& o, const aleph::containers::PointCloud<T>& pointCloud )
{
  auto n = pointCloud.size();
  for( decltype(n) i = 0; i < n; i++ )
  {
    auto data = pointCloud[i];

    for( auto it = data.begin(); it != data.end(); ++it )
    {
      if( it != data.begin() )
        o << "\t";

      o << *it;
    }

    o << "\n";
  }

  return o;
}

#endif
