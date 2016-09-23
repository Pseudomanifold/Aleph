#ifndef ALEPH_DISTANCES_DETAIL_MUNKRES_HH__
#define ALEPH_DISTANCES_DETAIL_MUNKRES_HH__

#include "Matrix.hh"

#include <algorithm>
#include <limits>
#include <vector>

namespace aleph
{

namespace distances
{

namespace detail
{

template <class T> class Munkres
{
public:

  void solve( Matrix<T>& matrix )
  {
    auto n = matrix.n();

    // Copy input matrix
    this->matrix = m;

    _rowMask = std::vector<bool>( size, false );
    _colMask = std::vector<bool>( size, false );

    subtractRowMinimum( matrix );

    unsigned short step = 1;

    while( step )
    {
      std::size_t row = 0;
      std::size_t col = 0;

      switch ( step ) {
        case 1:
          step = step1( matrix );
          break;
        case 2:
          step = step2( matrix );
          break;
        case 3:
          step = step3( matrix, row, col );
          break;
        case 4:
          step = step4( matrix );
          break;
        case 5:
          step = step5( matrix );
          break;
      }
    }

    for( std::size_t row = 0; row < n; row++ )
    {
      for( std::size_t col = 0; col < n; col++ )
      {
        if( _stars( row, col ) )
          matrix( row, col ) = T( 0 );
        else
          matrix( row, col ) = std::numeric_limits<T>::max();
      }
    }
  }

private:

  // Modifies matrix in-place by subtracting the minimum value in each
  // row. The function assumes that matrix is not empty.
  void subtractRowMinimum( Matrix<T>& matrix )
  {
    auto n = matrix.n();

    for( std::size_t row = 0; row < n; row++ )
    {
      auto min = matrix( row, 0 );

      for( std::size_t col = 0; col < n; col++ )
        min = std::min( min, matrix( row, col ) );

      for( std::size_t col = 0; col < n; col++ )
        matrix( row, col ) -= min;
    }
  }

  bool findUncoveredZeroInMatrix( const Matrix<T>& matrix, std::size_t& row, std::size_t& column ) const
  {
    auto n = matrix.n();

    for( row = 0; row < n; row++ )
    {
      if( !_rowMask[row] )
      {
        for( col = 0; col < n; col++ )
        {
          if( !_colMask[col] )
          {
            if( matrix( row, col ) == T( 0 ) )
              return true;
          }
        }
      }
    }

    return false;
  }

  // Find zeroes in the current matrix. If there is no starred zero in
  // the row or column, star the current value.
  unsigned short step1( Matrix<T>& matrix )
  {
    auto n = matrix.n();

    for( std::size_t row = 0; row < n; row++ )
    {
      for( std::size_t col = 0; col < n; col++ )
      {
        if( matrix( row, col ) == T( 0 ) )
        {
          // Check whether another zero in the same _column_ is already
          // starred.
          for( std::size_t r = 0; r < n; r++ )
          {
            if( _stars( r, col ) )
              goto skipCurrentColumn;
          }

          // Check whether another zero in the same _row_ is already
          // starred.
          for( std::size_t c = 0; c < n; c++ )
          {
            if( _stars( row, c ) )
              goto skipCurrentRow;
          }

          _stars( row, col )  = true;
          _primes( row, col ) = false;
        }
        skipCurrentColumn:
          ;
      }

      skipCurrentRow:
        ;
    }

    return 2;
  }

  // Cover each column that contains a starred zero. If enough columns
  // have been covered, the starred zeroes give us the complete set of
  // assignments.
  unsigned short step2( Matrix<T>& matrix )
  {
    auto n                     = matrix.n();
    std::size_t coveredColumns = 0;

    for( std::size_t row = 0; row < n; row++ )
    {
      for( std::size_t col = 0; col < n; col++ )
      {
        if( _stars( row, col ) )
        {
          _colMask[col] = true;
          ++coveredColumns;
        }
      }
    }

    if( coveredColumns >= n )
      return 0;

    return 3;
  }

  // Find an uncovered zero and prime it. If there is no starred zero in
  // the row containing this primed zero, go to step 5. Otherwise, cover
  // this row and uncover the column that contains the starred zero.
  unsigned short step3( Matrix<T>& matrix, std::size_t& row, std::size_t& col )
  {
    if( findUncoveredZeroInMatrix( matrix, row, col ) )
    {
      _primes( row, col ) = true;
      _stars( row, col )  = false;
    }
    else
      return 5;

    auto n = matrix.n();

    for( std::size_t c = 0; c < n; c++ )
    {
      if( _stars( row, c ) )
      {
        _rowMask[row] = true;
        _colMask[col] = false;

        return 3;
      }
    }

    return 4;
  }

  unsigned short step4( Matrix<T>& matrix, std::size_t row, std::size_t col )
  {
    auto n = matrix.n();

    using Pair     = std::pair<std::size_t, std::size_t>;
    using Sequence = std::vector<Pair>;

    Sequence sequence;

    sequence.push_back( std::make_pair( row, col ) );

    Pair p1 = std::make_pair( 0, 0 );
    Pair p2 = std::make_pair( 0, 0 );

    std::size_t r = 0;
    std::size_t c = col;

    size_t row, col = savecol;

    bool havePair = false;
    do
    {
      havePair = false;
      for ( r = 0; r < n; r++ )
      {
        if( _stars(r,c) )
        {
          p1.first  = r;
          p1.second = c;

          if( std::find( sequence.begin(), sequence.end(), p1 ) != sequence.end() )
            continue;

          havePair = true;

          sequence.push_back( p1 );
          break;
        }
      }

      if( !havePair )
        break;

      havePair = false;

      for( c = 0; c < n; c++ )
      {
        if( _primes(r,c) )
        {
          p2.first  = r;
          p2.second = c;

          if( std::find( sequence.begin(), sequence.end(), p2 ) != sequence.end() )
            continue;

          havePair = true;

          sequence.push_back( p2 );
          break;
        }
      }
    }
    while ( havePair );

    for( auto&& pair : sequence )
    {
      // Un-star
      if( _stars(pair.first, pair.second) )
        _stars(pair.first, pair.second) = false;

      // Star each primed zero
      if( _primes(pair.first, pair.second) )
        _stars(pair.first, pair.second) = true;
    }

    // Erase all primes & uncover all columns and rows -----------------

    for( std::size_t row = 0; row < n; row++ )
    {
      for( std::size_t col = 0; col < n; col++ )
        _primes(row, col) = false;

      _rowMask[row] = false;
      _colMask[row] = false;
    }

    return 2;
  }

  unsigned short step5( Matrix<T>& matrix )
  {
    auto n = matrix.n();
    T v    = std::numeric_limits<T>::max();

    for( std::size_t row = 0; row < n; row++ )
    {
      if( !_rowMask[row] )
      {
        for( std::size_t col = 0; col < n; col++ )
        {
          if( !_colMask[col] )
            if( matrix( row, col ) != T( 0 ) && matrix( row, col ) < v )
              v = matrix( row, col );
        }
      }
    }

    for( std::size_t row = 0; row < n; row++ )
    {
      if( _rowMask[row] )
      {
        for( std::size_t col = 0; col < n; col++ )
          matrix( row, col ) += v;
      }
    }

    for( std::size_t col = 0; col < n; col++ )
    {
      if( !_colMask[col] )
      {
        for( std::size_t row = 0; row < n; row++ )
          matrix( row, col ) -= v;
      }
    }

    return 3;
  }

  Matrix<bool> _stars;
  Matrix<bool> _primes;

  std::vector<bool> _rowMask;
  std::vector<bool> _colMask;
};


}

}

}

#endif
