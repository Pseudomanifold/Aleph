#ifndef ALEPH_DISTANCES_DETAIL_MATRIX_HH__
#define ALEPH_DISTANCES_DETAIL_MATRIX_HH__

#include <algorithm>

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

  Matrix( const Matrix& other )           = delete;
  Matrix operator=( const Matrix& other ) = delete;

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

}

}

}

#endif
