#ifndef ALEPH_DISTANCES_DETAIL_MUNKRES_HH__
#define ALEPH_DISTANCES_DETAIL_MUNKRES_HH__

#include <aleph/distances/detail/Matrix.hh>

#include <algorithm>
#include <limits>
#include <numeric>
#include <utility>
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

  /**
    Given a reduced matrix stored by the class and a matrix of costs,
    calculates the costs and returns them. Note that this function is
    unable to know whether the stored matrix has been reduced or not.
  */

  T cost( const Matrix<T>& costs ) const noexcept
  {
    std::vector<T> allCosts;

    auto n = _matrix.n();

    if( n != costs.n() )
    {
      if( std::numeric_limits<T>::has_quiet_NaN )
        return std::numeric_limits<T>::quiet_NaN();
      else
        return std::numeric_limits<T>::max();
    }

    for( std::size_t row = 0; row < n; row++ )
    {
      for( std::size_t col = 0; col < n; col++ )
      {
        // This value denotes that a match has been made between the
        // current row and the current column.
        if( _matrix( row, col ) == T() )
          allCosts.push_back( costs( row, col ) );
      }
    }

    return std::accumulate( allCosts.begin(), allCosts.end(), T() );
  }

  /**
    Stores the matching in an output iterator. The output iterator needs
    to be able to handle pairs of indices. Note that this function needs
    an
  */

  template <class OutputIterator> void matching( OutputIterator result ) const noexcept
  {
    using Index = std::size_t;
    auto n      = _matrix.n();

    for( std::size_t row = 0; row < n; row++ )
    {
      for( std::size_t col = 0; col < n; col++ )
      {
        if( _matrix( row, col ) == T() )
          *result++ = std::make_pair( static_cast<Index>( row ), static_cast<Index>( col ) );
      }
    }
  }

  // Constructor -------------------------------------------------------

  Munkres( const Matrix<T>& matrix )
    : _matrix( matrix )
    , _stars( matrix.n() )
    , _primes( matrix.n() )
    , _rowMask( std::vector<bool>( _matrix.n(), false ) )
    , _colMask( std::vector<bool>( _matrix.n(), false ) )
  {
  }

  // Solver ------------------------------------------------------------

  Matrix<T> operator()()
  {
    subtractRowMinimum( _matrix );

    unsigned short step = 1;

    std::size_t row = 0;
    std::size_t col = 0;

    while( step != 0 )
    {
      switch( step )
      {
      case 1:
        step = step1();
        break;
      case 2:
        step = step2();
        break;
      case 3:
        step = step3( _matrix, row, col );
        break;
      case 4:
        step = step4( _matrix, row, col );
        break;
      case 5:
        step = step5( _matrix );
        break;
      }
    }

    // Prepare reduced matrix ------------------------------------------

    auto n = _matrix.n();

    for( std::size_t row = 0; row < n; row++ )
    {
      for( std::size_t col = 0; col < n; col++ )
      {
        if( _stars( row, col ) )
          _matrix( row, col ) = T( 0 );
        else
          _matrix( row, col ) = std::numeric_limits<T>::max();
      }
    }

    return _matrix;
  }

private:

  // Modifies matrix in-place by subtracting the minimum value in each
  // row. The function assumes that matrix is not empty.
  static void subtractRowMinimum( Matrix<T>& matrix ) noexcept
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

  bool findUncoveredZeroInMatrix( std::size_t& row, std::size_t& col ) const
  {
    auto n = _matrix.n();

    for( row = 0; row < n; row++ )
    {
      if( !_rowMask[row] )
      {
        for( col = 0; col < n; col++ )
        {
          if( !_colMask[col] )
          {
            if( _matrix( row, col ) == T( 0 ) )
              return true;
          }
        }
      }
    }

    return false;
  }

  // Find zeroes in the current matrix. If there is no starred zero in
  // the row or column, star the current value.
  unsigned short step1()
  {
    auto n = _matrix.n();

    for( std::size_t row = 0; row < n; row++ )
    {
      for( std::size_t col = 0; col < n; col++ )
      {
        if( _matrix( row, col ) == T( 0 ) )
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
  unsigned short step2()
  {
    auto n                     = _matrix.n();
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
    if( findUncoveredZeroInMatrix( row, col ) )
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
        _colMask[c]   = false;

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
      {
        _stars(pair.first, pair.second)  = false;
        _primes(pair.first, pair.second) = false;
      }

      // Star each primed zero
      if( _primes(pair.first, pair.second) )
      {
        _stars(pair.first, pair.second)  = true;
        _primes(pair.first, pair.second) = false;
      }
    }

    // Erase all primes & uncover all columns and rows -----------------

    for( std::size_t row = 0; row < n; row++ )
    {
      for( std::size_t col = 0; col < n; col++ )
      {
        if( _primes( row, col ) )
        {
          _primes(row, col) = false;
          _stars(row, col)  = false;
        }
      }

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

  // Attributes --------------------------------------------------------

  // The matrix is going to be reduced in-place. The class does not
  // store additional information about its status. Hence, there is
  // no way of determining whether the algorithm has converged.
  Matrix<T> _matrix;

  Matrix<bool> _stars;
  Matrix<bool> _primes;

  std::vector<bool> _rowMask;
  std::vector<bool> _colMask;
};

}

}

}

#endif
