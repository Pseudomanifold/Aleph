#ifndef ALEPH_MATH_SYMMETRIC_MATRIX_HH__
#define ALEPH_MATH_SYMMETRIC_MATRIX_HH__

#include <algorithm>
#include <iomanip>
#include <iosfwd>
#include <stdexcept>

#include <cstddef>

namespace aleph
{

namespace math
{

/**
  @class SymmetricMatrix

  Describes a symmetric matrix of an arbitrary data type. The class is
  able to access to the data using an interface that pretends to model
  a regular matrix.

  Hence both types of queries will be valid:

      M(i,j)
      M(j,i)

  Other than that, this class aims to have a small footprint.

  @tparam T Data type stored in matrix, e.g. `double`
  @tparam I Index type for accessing the matrix. You may change this for
            small matrices in order to reduce their memory usage. Notice
            that the default values should be sufficient in most cases.
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
    Creates a symmetric matrix as a copy of another symmetric matrix, by
    simply copying all data values. This function first allocates memory
    before clearing it, and finally fills it with a copy of the data the
    other matrix stores. Hence, this function is *not* highly-efficient.

    @param other Matrix from which to copy data
  */

  SymmetricMatrix( const SymmetricMatrix& other )
    : SymmetricMatrix( other._numRows )
  {
    std::copy( other._data, other._data + other._size, _data );
  }

  /** Destroys the symmetric matrix */
  ~SymmetricMatrix()
  {
    delete[] _data;
  }

  /** Swaps two matrices */
  void swap( SymmetricMatrix& other )
  {
    std::swap( _numRows, other._numRows );
    std::swap( _size   , other._size    );
    std::swap( _data   , other._data    );
  }

  /**
    Assigns a given symmetric matrix to the current symmetric matrix,
    taking care not to throw an error.
  */

  SymmetricMatrix& operator=( SymmetricMatrix other )
  {
    this->swap( other );
    return *this;
  }

  /**
    Provides element-wise access to the matrix and returns the element
    at the specified position. The function throws if an invalid index
    is encountered.

    @param row    Desired row for query
    @param column Desired column for query

    @returns Const reference to the desired value
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

  /** Returns number of rows */
  I numRows() const noexcept
  {
    return _numRows;
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

/**
  Output operator for pretty-printing symmetric matrices. This is most
  useful for debugging. The function takes care of formatting but will
  not modify the flags of the output stream permanently.

  @param o Output stream
  @param M Input matrix

  @returns Modified output stream
*/

template <class T> std::ostream& operator<<( std::ostream& o, const aleph::math::SymmetricMatrix<T>& M )
{
  auto flags( o.flags() );

  auto n = M.numRows();

  o << std::setfill( ' ' );

  for( decltype(n) i = 0; i < n; i++ )
  {
    for( decltype(n) j = 0; j < n; j++ )
    {
      if( j != 0 )
        o << ", ";

      o << std::setw( 10 ) << std::right << M(i,j);
    }

    o << "\n";
  }

  o.flags( flags );
  return o;
}

#endif
