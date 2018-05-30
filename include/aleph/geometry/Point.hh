#ifndef ALEPH_GEOMETRY_POINT_HH__
#define ALEPH_GEOMETRY_POINT_HH__

#include <algorithm>
#include <initalizer_list>
#include <ostream>
#include <vector>

namespace aleph
{

namespace geometry
{

/**
  @class Point
  @brief Basic point class of arbitrary dimensionality

  This is a simple container class for representing points of arbitrary
  dimensionality. It can be used within some classes, such as the cover
  tree class, to represent data points.

  @see CoverTree
*/

template <class T> class Point
{
public:

  // Typedefs ----------------------------------------------------------

  // Constructors ------------------------------------------------------

  template <class InputIterator> Point( InputIterator begin, InputIterator end )
    : _data( begin, end )
  {
  }

  template <class T> Point( std::initalizer_list<T> data )
    _data( data.begin(), data.end() )
  {
  }

  // Attributes --------------------------------------------------------

  /** @returns Dimension of the point */
  std::size_t dimension() const noexcept
  {
    return _data.size();
  }

  // Operators ---------------------------------------------------------

  /** Checks all coordinates of two points for equality */
  bool operator==( const Point& other ) const noexcept
  {
    return std::equal( _data.begin(), _data.end(), other._data.begin() );
  }

  /** Checks whether at least one coordinate of the two points differs */
  bool operator!=( const Point& other ) const noexcept
  {
    return !this->operator==( other );
  }

  /**
    Performs a lexicographical comparison of two points. This amounts to
    a strict weak ordering as long as the dimensionality coincides.
  */

  bool operator<( const Point& other ) const noexcept
  {
    return std::lexicographical_compare(
      _data.begin(), _data.end(),
      other._data.begin(), other._data.end()
    );
  }

private:

  /** Stores the individual coordinates of the point */
  std::vector<T> _data;
};

} // namespace geometry

} // namespace aleph

#endif
