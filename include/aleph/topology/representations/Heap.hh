#ifndef ALEPH_REPRESENTATIONS_HEAP_HH__
#define ALEPH_REPRESENTATIONS_HEAP_HH__

#include <algorithm>
#include <utility>
#include <vector>

namespace aleph
{

namespace topology
{

namespace representations
{

template <class IndexType = unsigned> class Heap
{
public:
  using Index = IndexType;

  void setNumColumns( Index numColumns )
  {
    _data.resize( static_cast<std::size_t>( numColumns ) );
    _dimensions.resize( static_cast<std::size_t>( numColumns ) );
  }

  Index getNumColumns() const
  {
    return static_cast<Index>( _data.size() );
  }

  std::pair<Index, bool> getMaximumIndex( Index column ) const
  {
    if( _data.at( static_cast<std::size_t>( column ) ).empty() )
      return std::make_pair( Index(0), false );
    else
    {
      auto&& column_ = _data.at( static_cast<std::size_t>( column ) );
      auto index     = column_.front();

      std::pop_heap( column_.begin(), column_.end() );

      while( !column_.empty() && column_.front() == index )
      {
        std::pop_heap( column_.begin(), column_.end() );
        column_.pop_back();

        if( column_.empty() )
          return std::make_pair( Index(0), false );
      }

      return std::make_pair( index, true );
    }
  }

  void addColumns( Index source, Index target )
  {
    auto&& sourceColumn = _data.at( static_cast<std::size_t>( source ) );
    auto&& targetColumn = _data.at( static_cast<std::size_t>( target ) );

    targetColumn.reserve( sourceColumn.size() + targetColumn.size() );

    for( auto&& value : sourceColumn )
    {
      targetColumn.push_back( value );
      std::push_heap( targetColumn.begin(), targetColumn.end() );
    }
  }

  template <class InputIterator> void setColumn( Index column,
                                                 InputIterator begin, InputIterator end )
  {
    _data.at( static_cast<std::size_t>( column ) ).assign( begin, end );

    // Ensures proper heap order. Else, the reduction algorithm will
    // not be able to reduce the matrix.
    std::make_heap( _data.at( static_cast<std::size_t>( column ) ).begin(), _data.at( static_cast<std::size_t>( column ) ).end() );

    // Upon initialization, the column must by necessity have the dimension
    // that is indicated by the amount of indices in its boundary. The case
    // of 0-simplices needs special handling.
    _dimensions.at( static_cast<std::size_t>( column ) )
        = begin == end ? 0
                       : static_cast<Index>( std::distance( begin, end ) - 1 );
  }

  std::vector<Index> getColumn( Index column ) const
  {
    return _data.at( static_cast<std::size_t>( column ) );
  }

  void clearColumn( Index column )
  {
    _data.at( static_cast<std::size_t>( column ) ).clear();
  }

  void setDimension( Index column, Index dimension )
  {
    _dimensions.at( static_cast<std::size_t>( column ) ) = dimension;
  }

  Index getDimension( Index column ) const
  {
    return _dimensions.at( static_cast<std::size_t>( column ) );
  }

  Index getDimension() const
  {
    if( _dimensions.empty() )
      return Index(0);
    else
      return *std::max_element( _dimensions.begin(), _dimensions.end() );
  }

  bool operator==( const Heap& other ) const
  {
    return _data == other._data && _dimensions == other._dimensions;
  }

private:
  std::vector< std::vector<Index> > _data;
  std::vector<Index> _dimensions;
};

} // namespace representations

} // namespace topology

} // namespace aleph

#endif
