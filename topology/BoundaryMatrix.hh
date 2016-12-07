#ifndef ALEPH_BOUNDARY_MATRIX_HH__
#define ALEPH_BOUNDARY_MATRIX_HH__

#include <algorithm>
#include <fstream>
#include <istream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace aleph
{

namespace topology
{

template <class Representation> class BoundaryMatrix
{
public:
  using Index = typename Representation::Index;

  void setNumColumns( Index numColumns )
  {
    _representation.setNumColumns( numColumns );
  }

  Index getNumColumns() const
  {
    return _representation.getNumColumns();
  }

  std::pair<Index, bool> getMaximumIndex( Index column ) const
  {
    return _representation.getMaximumIndex( column );
  };

  void addColumns( Index source, Index target )
  {
    _representation.addColumns( source, target );
  };

  template <class InputIterator> void setColumn( Index column,
                                                 InputIterator begin, InputIterator end )

  {
    _representation.setColumn( column, begin, end );
  }

  std::vector<Index> getColumn( Index column ) const
  {
    return _representation.getColumn( column );
  }

  void clearColumn( Index column )
  {
    _representation.clearColumn( column );
  }

  void setDimension( Index column, Index dimension )
  {
    _representation.setDimension( column, dimension );
  }

  Index getDimension( Index column ) const
  {
    return _representation.getDimension( column );
  }

  Index getDimension() const
  {
    return _representation.getDimension();
  }

  bool isDualized() const
  {
    return _isDualized;
  }

  // Comparison --------------------------------------------------------

  bool operator==( const BoundaryMatrix& other ) const
  {
    return    _representation == other._representation
           && _isDualized     == other._isDualized;
  }

  bool operator!=( const BoundaryMatrix& other ) const
  {
    return !this->operator==( other );
  }

  // Dualization -------------------------------------------------------

  BoundaryMatrix dualize() const
  {
    auto&& numColumns = this->getNumColumns();

    std::vector< std::vector<Index> > dualMatrix( numColumns );
    std::vector<Index> dualDimensions( numColumns );
    std::vector<std::size_t> dualColumnSizes( numColumns );

    // Determine the size of every column in the dualized matrix. This
    // keeps memory re-allocation at a minimum.

    for( Index j = 0; j < numColumns; j++ )
    {
      auto&& column = this->getColumn(j);

      for( auto&& i : column )
        ++dualColumnSizes[ numColumns - 1 - i ];
    }

    for( Index j = 0; j < numColumns; j++ )
      dualMatrix[j].reserve( dualColumnSizes[j] );

    // Calculate the actual anti-transpose of the matrix. Since the
    // vectors have been properly resized, this operation should be
    // relatively well-behaved.

    for( Index j = 0; j < numColumns; j++ )
    {
      auto&& column = this->getColumn( j );

      for( auto&& i : column )
        dualMatrix[ numColumns - 1 - i ].push_back( numColumns - 1 - j );
    }

    auto&& d = this->getDimension();

    // FIXME: Do I need this?
    for( Index j = 0; j < numColumns; j++ )
      dualDimensions.at( numColumns - 1 - j ) = d - this->getDimension( j ); // FIXME: Change operator later after debugging

    BoundaryMatrix<Representation> M;
    M.setNumColumns( static_cast<Index>( dualMatrix.size() ) );

    for( Index j = 0; j < M.getNumColumns(); j++ )
    {
      // Do not assume that the column is properly sorted. A normal
      // std::reverse should be sufficient in most cases here but I
      // do not want to take any chances.
      std::sort( dualMatrix.at(j).begin(), dualMatrix.at(j).end() );

      M.setColumn( j,
                   dualMatrix.at(j).begin(), dualMatrix.at(j).end() );

      M.setDimension( j, dualDimensions.at(j) );
    }

    M._isDualized = !this->isDualized();
    return M;
  }

  // I/O operations ----------------------------------------------------

  static BoundaryMatrix load( std::istream& in )
  {
    std::string line;
    Index numColumns = Index(0);

    while( std::getline( in, line ) )
    {
      // Remove all spaces
      line.erase( std::remove_if( line.begin(), line.end(), ::isspace ),
                  line.end() );

      if( line.empty() == false && line.front() != '#' )
        ++numColumns;
    }

    in.clear();
    in.seekg( 0 );

    BoundaryMatrix<Representation> M;
    M.setNumColumns( numColumns );

    Index curColumn = Index(0);

    while( std::getline( in, line ) )
    {
      // Ignore empty lines and comment lines
      if( line.empty() || line.front() == '#' )
        continue;

      std::istringstream iss( line );

      std::vector<std::string> tokens;

      std::copy( std::istream_iterator<std::string>( iss ),
                 std::istream_iterator<std::string>(),
                 std::back_inserter( tokens ) );

      // Ignore this line if it does not contain any tokens
      if( tokens.empty() )
        continue;

      std::vector<Index> indices;
      indices.reserve( tokens.size() );

      std::transform( tokens.begin(), tokens.end(), std::back_inserter( indices ),
                      [] ( const std::string& token )
                      {
                        return static_cast<Index>( std::stoul( token ) );
                      } );

      if( indices.empty() )
        throw std::runtime_error( "Amount of indices in boundary must not be empty" );

      // TODO: This ignores the dimension of the column
      // TODO: This assumes that the column indices are ordered
      M.setColumn( curColumn,
                   indices.begin() + 1, indices.end() );

      if( M.getDimension( curColumn ) != indices.front() )
        throw std::runtime_error( "Inconsistency between actual number of indices and specified number of indices in boundary" );

      ++curColumn;
    }

    return M;
  }

  static BoundaryMatrix load( const std::string& filename )
  {
    std::ifstream in( filename );
    return BoundaryMatrix::load( in );
  }

private:
  Representation _representation;

  /**
    Flag indicating whether the matrix is dualized or not. By default
    no matrix is dualized. This flag is used by some of the reduction
    algorithms to determine how to calculate indices.
  */

  bool _isDualized = false;
};

// ---------------------------------------------------------------------

template <class Representation> std::ostream& operator<< ( std::ostream& o, const BoundaryMatrix<Representation>& M )
{
  using Index = typename Representation::Index;

  auto numColumns = M.getNumColumns();

  for( Index j = Index(0); j < numColumns; ++j )
  {
    auto column = M.getColumn( j );

    if( !column.empty() )
    {
      for( auto&& c : column )
        o << c << " ";
    }
    else
      o << "-";

    o << "\n";
  }

  return o;
}

// ---------------------------------------------------------------------

} // namespace topology

} // namespace aleph

#endif
