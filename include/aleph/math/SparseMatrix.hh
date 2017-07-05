#ifndef ALEPH_MATH_SPARSE_MATRIX_HH__
#define ALEPH_MATH_SPARSE_MATRIX_HH__

#include <cstddef>

#include <algorithm>
#include <set>
#include <stdexcept>
#include <vector>

namespace aleph
{

namespace math
{

template <class T> class SparseBinaryMatrix
{
public:
  using IndexType   = T;
  using ColumnType  = std::set<IndexType>;
  using ColumnsType = std::vector<ColumnType>;
  using IndexMap    = std::vector<IndexType>;
  using size_type   = typename ColumnType::size_type;

  explicit SparseBinaryMatrix( IndexType columns )
  {
    _columns.resize( size_type( columns ) );
  }

  /**
    Sets the value of a given entry in the matrix. If the client
    specifies that the value of the entry shall be true, it will
    be added to the column. Else, the entry will be removed.
  */

  void set( IndexType row, IndexType column )
  {
    _columns.at( size_type( column ) ).insert( row );
  }

  /**
    Gets the value of a given entry in the matrix. Calling this
    function is well-defined, even for non-existent entries. It
    will throw, however, when an invalid column index is used.
  */

  bool get( IndexType row, IndexType column ) const
  {
    auto&& c = _columns.at( size_type( column ) );
    return c.find( row ) != c.end();
  }

  /** Returns all non-zero values in a given column. */
  template <class OutputIterator> void get( IndexType column,
                                            OutputIterator result ) const
  {
    auto&& c = _columns.at( size_type( column ) );

    std::copy( c.begin(), c.end(), result );
  }

  /** Sets row/column indices */
  template <class InputIterator> void setIndices( InputIterator begin, InputIterator end )
  {
    _indexMap.assign( begin, end );
    if( _indexMap.size() != this->numColumns() )
      throw std::length_error( "Number of indices does not coincide with number of columns" );
  }

  /**
    Gets row/class index of a given column. If no foreign indices have
    been set by the client, this function just returns its input.
  */

  IndexType getIndex( IndexType column ) const
  {
    if( _indexMap.empty() )
      return column;
    else
      return _indexMap.at( size_type( column ) );
  }

  /** Returns number of non-zero entries in a given column */
  std::size_t numEntries( IndexType column ) const
  {
    return _columns.at( size_type( column ) ).size();
  }

  /** Returns number of columns */
  std::size_t numColumns() const noexcept
  {
    return _columns.size();
  }

private:

  // The index map maps a local row/column index in [0..n-1] to another
  // index that has been specified by the client. This is used whenever
  // indices occur that may not be contiguous.
  IndexMap _indexMap;
  ColumnsType _columns;
};

template <class I, class D> class SparseMatrix
{
  /*
    TODO: Implement this. In contrast to the binary matrix, this matrix
    should also be capable of storing the value in a given column.

    Hence, there is the need for an additional lookup data structure.
  */
};

} // namespace math

} // namespace aleph

#endif
