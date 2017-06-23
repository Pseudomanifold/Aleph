#ifndef ALEPH_DISTANCES_DETAIL_MATRIX_HH__
#define ALEPH_DISTANCES_DETAIL_MATRIX_HH__

#include <algorithm>
#include <ostream>
#include <utility>

#include <cassert>
#include <cstddef>

namespace aleph
{

namespace distances
{

namespace detail
{

template <class T> class Matrix
{
public:

  Matrix( std::size_t n )
    : _n( n )
    , _data( new T*[_n] )
  {
    for( std::size_t row = 0; row < _n; row ++ )
    {
      _data[row] = new T[_n];

      std::fill( _data[row], _data[row] + _n, T() );
    }
  }

  ~Matrix()
  {
    for( std::size_t row = 0; row < _n; row++ )
      delete[] _data[row];

    delete[] _data;
  }

  Matrix( const Matrix& other )
    : _n( other._n )
    , _data( new T*[_n] )
  {
    for( std::size_t row = 0; row < _n; row ++ )
    {
      _data[row] = new T[_n];

      std::copy( other._data[row], other._data[row] + _n, _data[row] );
    }
  }

  Matrix( Matrix&& other )
    : Matrix( 0 )
  {
    swap( *this, other );
  }

  Matrix& operator=( Matrix other )
  {
    swap( *this, other );
    return *this;
  }

  friend void swap( Matrix& m1, Matrix& m2 ) noexcept
  {
    using std::swap;

    swap( m1._n,    m2._n );
    swap( m1._data, m2._data );
  }

  std::size_t n() const
  {
    return _n;
  }

  T& operator()( std::size_t row, std::size_t column )
  {
    return const_cast<T&>( static_cast<const Matrix&>( *this ).operator()( row, column ) );
  }

  const T& operator()( std::size_t row, std::size_t column ) const
  {
    assert( row    < _n );
    assert( column < _n );

    return _data[row][column];
  }

private:

  std::size_t _n = 0;
  T** _data      = nullptr;
};

template <class T> std::ostream& operator<<( std::ostream& o, const Matrix<T>& m )
{
  auto n = m.n();

  for( decltype(n) row = 0; row < n; row++ )
  {
    for( decltype(n) col = 0; col < n; col++ )
    {
      if( col > 0 )
        o << " ";

      o << m(row,col);
    }

    o << "\n";
  }

  return o;
}

}

}

}

#endif
