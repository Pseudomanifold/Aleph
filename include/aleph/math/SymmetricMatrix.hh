#ifndef ALEPH_MATH_SYMMETRIC_MATRIX_HH__
#define ALEPH_MATH_SYMMETRIC_MATRIX_HH__

#include <algorithm>
#include <stdexcept>

#include <cstddef>

namespace aleph
{

namespace math
{

/**
  :class SymmetricMatrix: Describes a symmetric matrix of an arbitrary
  data type. This class provides access to the data using an interface
  that pretends to model a regular matrix.

  Hence both types of queries will be valid:::

    M(i,j)
    M(j,i)

  Other than that, this class aims to have a small footprint.
*/

template <class T, class I = std::size_t> class SymmetricMatrix
{
public:

  /** Creates an empty symmetric matrix */
  SymmetricMatrix()
  {
  }

  /** Creates symmetric matrix with a given number of rows and columns */
  SymmetricMatrix( I n )
    : _numRows( n )
    , _size( n * ( n + 1 ) / 2 )
    , _data( new T[ _size ] )
  {
    std::fill( _data, _data + _size, I(0) );
  }

  /**
    Provides element-wise access to the matrix and returns the element
    at the specified position. The function throws if an invalid index
    is encountered.
  */

  const T& operator()( I row, I column ) const
  {
    if( row >= _numRows || column >= _numRows )
      throw std::out_of_range( "Index is out of range" );

    if( row > column )
      std::swap( row, column );

    auto index = row * _numRows - row * ( row + 1 ) / 2 + column;

    if( index >= _size )
      throw std::out_of_range( "Index is out of range" );
    else
      return _data[index];
  }

  T& operator()( I row, I column )
  {
    return const_cast<T&>( static_cast<const SymmetricMatrix&>( *this )( row, column ) );
  }

  /** Returns size of matrix */
  I size() const noexcept
  {
    return _size;
  }

  /** Checks whether the matrix is empty */
  bool empty() const noexcept
  {
    return _numRows == 0 && _size == 0 && _data == nullptr;
  }

private:

  /** Number of rows (and columns) */
  I _numRows = I(0);

  /** Total size of data storage */
  I _size = I(0);

  /** 1D data storage (for efficiency reasons) */
  T* _data = nullptr;

};

} // namespace math

} // namespace aleph

#endif