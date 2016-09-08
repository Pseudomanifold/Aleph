#ifndef ALEPH_BOUNDARY_MATRIX_HH__
#define ALEPH_BOUNDARY_MATRIX_HH__

#include <vector>

namespace aleph
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

private:
  Representation _representation;
};

}

#endif
