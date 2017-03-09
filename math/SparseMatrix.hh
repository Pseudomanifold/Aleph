#ifndef ALEPH_MATH_SPARSE_MATRIX_HH__
#define ALEPH_MATH_SPARSE_MATRIX_HH__

#include <cstddef>

#include <algorithm>
#include <set>
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


  SparseBinaryMatrix( IndexType columns )
  {
    _columns.resize( columns );
  }

  /**
    Sets the value of a given entry in the matrix. If the client
    specifies that the value of the entry shall be true, it will
    be added to the column. Else, the entry will be removed.
  */

  void set( IndeType row, IndexType column )
  {
    _columns.at( column ).insert( row );
  }

  /**
    Gets the value of a given entry in the matrix. Calling this
    function is well-defined, even for non-existent entries. It
    will throw, however, when an invalid column index is used.
  */

  bool get( IndexType row, IndexType column ) const
  {
    auto&& column = _columns.at( column );
    return column.find( row ) != column.end();
  }

  /** Returns all non-zero values in a given column. */
  template <class OutputIterator> void get( IndexType column,
                                            OutputIterator result ) const
  {
    auto&& column = _columns.at( column );

    std::copy( column.begin(), column.end(), result );
  }

  std::size_t numColumns() const noexcept
  {
    return _columns.size();
  }


private:

  ColumnsType _columns;
};

} // namespace math

} // namespace aleph

#endif
