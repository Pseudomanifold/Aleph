#ifndef ALEPH_PERSISTENCE_DIAGRAM_HH__
#define ALEPH_PERSISTENCE_DIAGRAM_HH__

#include <limits>
#include <utility>
#include <vector>

namespace aleph
{

template <class DataType> class PersistenceDiagram
{
public:

  class Point
  {
  public:

    Point( DataType x )
      : _x( x )
      , _y( std::numeric_limits<DataType>::max() )
    {
      if( std::numeric_limits<DataType>::has_infinity )
        _y = std::numeric_limits<DataType>::infinity();
    }

    Point( DataType x, DataType y )
      : _x( x )
      , _y( y )
    {
    }

    DataType x() const { return _x; }
    DataType y() const { return _y; }

    DataType persistence() const
    {
      return _y - _x;
    }

  private:
    DataType _x;
    DataType _y;
  };

  // Typedefs & aliases ------------------------------------------------

  using ValueType     = Point;
  using ContainerType = std::vector<ValueType>;

  using ConstIterator = typename ContainerType::const_iterator;
  using Iterator      = typename ContainerType::iterator;

  using value_type    = ValueType;

  // Iterators ---------------------------------------------------------

  ConstIterator begin() const { return _points.begin(); }
  Iterator      begin()       { return _points.begin(); }

  ConstIterator end() const   { return _points.end(); }
  Iterator      end()         { return _points.end(); }

  // Modification ------------------------------------------------------

  void add( DataType x )
  {
    _points.push_back( Point( x ) );
  }

  void add( DataType x, DataType y )
  {
    _points.push_back( Point( x, y ) );
  }

  Iterator erase( Iterator position )
  {
    return _points.erase( position );
  }

  Iterator erase( Iterator begin, Iterator end )
  {
    return _points.erase( begin, end );
  }

  // Queries -----------------------------------------------------------

  std::size_t dimension() const
  {
    return _dimension;
  }

  std::size_t size() const
  {
    return _points.size();
  }

  bool empty() const
  {
    return _points.empty();
  }

private:

  /** Dimension of the persistence pairs stored in the diagram */
  std::size_t _dimension = 0;

  /** Container of persistence pairs */
  ContainerType _points;
};

}

#endif
