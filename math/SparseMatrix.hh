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

  void set( IndexType row, IndexType column )
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
    auto&& c = _columns.at( column );
    return c.find( row ) != c.end();
  }

  /** Returns all non-zero values in a given column. */
  template <class OutputIterator> void get( IndexType column,
                                            OutputIterator result ) const
  {
    auto&& c = _columns.at( column );

    std::copy( c.begin(), c.end(), result );
  }

  std::size_t numColumns() const noexcept
  {
    return _columns.size();
  }


private:

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
