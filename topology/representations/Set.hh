#ifndef ALEPH_REPRESENTATIONS_SET_HH__
#define ALEPH_REPRESENTATIONS_SET_HH__

#include <algorithm>
#include <set>
#include <vector>

namespace aleph
{

namespace topology
{

namespace representations
{

template <class IndexType> class Set
{
public:
  using Index = IndexType;

  void setNumColumns( Index numColumns )
  {
    _data.resize( numColumns );
    _dimensions.resize( numColumns );
  }

  Index getNumColumns() const
  {
    return static_cast<Index>( _data.size() );
  }

  std::pair<Index, bool> getMaximumIndex( Index column ) const
  {
    if( _data.at( column ).empty() )
      return std::make_pair( Index(0), false );
    else
      return std::make_pair( *_data.at( column ).rbegin(), true );
  }

  void addColumns( Index source, Index target )
  {
    auto&& sourceColumn = _data.at( source );
    auto&& targetColumn = _data.at( target );

    std::set<Index> result;

    std::set_symmetric_difference( sourceColumn.begin(), sourceColumn.end(),
                                   targetColumn.begin(), targetColumn.end(),
                                   std::inserter( result, result.begin() ) );

    targetColumn.swap( result );
  }

  template <class InputIterator> void setColumn( Index column,
                                                 InputIterator begin, InputIterator end )
  {
    _data.at( column ).clear();
    _data.at( column ).insert( begin, end );

    // Upon initialization, the column must by necessity have the dimension
    // that is indicated by the amount of indices in its boundary. The case
    // of 0-simplices needs special handling.
    _dimensions.at( column )
        = begin == end ? 0
                       : static_cast<Index>( std::distance( begin, end ) - 1 );
  }

  std::vector<Index> getColumn( Index column ) const
  {
    return { _data.at( column ).begin(), _data.at( column ).end() };
  }

  void clearColumn( Index column )
  {
    _data.at( static_cast<std::size_t>( column ) ).clear();
  }

  void setDimension( Index column, Index dimension )
  {
    _dimensions.at( column ) = dimension;
  }

  Index getDimension( Index column ) const
  {
    return _dimensions.at( column );
  }

  Index getDimension() const
  {
    if( _dimensions.empty() )
      return Index(0);
    else
      return *std::max_element( _dimensions.begin(), _dimensions.end() );
  }

  bool operator==( const Set& other ) const
  {
    return _data == other._data && _dimensions == other._dimensions;
  }

private:
  std::vector< std::set<Index> > _data;
  std::vector<Index> _dimensions;
};

} // namespace representations

} // namespace topology

} // namespace aleph

#endif
