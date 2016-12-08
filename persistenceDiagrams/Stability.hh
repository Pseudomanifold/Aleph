#ifndef ALEPH_PERSISTENCE_DIAGRAMS_STABILITY_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_STABILITY_HH__

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <algorithm>
#include <utility>
#include <vector>

// FIXME: Remove after debugging
#include <iostream>

namespace aleph
{

namespace detail
{

template <class T> struct EventPoint
{
  T value;
  bool destroyer;

  bool operator<( const EventPoint& other ) const noexcept
  {
    // Sort event point by their corresponding values. In case of ties,
    // consider creators to come after destroyers. This makes sense, as
    // a creator may increase the number of active intervals again.
    return value < other.value || ( value == other.value && destroyer && !other.destroyer );
  }
};

} // namespace detail

template <class DataType> std::vector< std::pair<DataType, DataType> > stabilityFunction( const PersistenceDiagram<DataType>& D )
{
  using namespace detail;
  using EP = EventPoint<DataType>;

  std::vector<EP> eventPoints;
  eventPoints.reserve( 2 * D.size() );

  for( auto&& p : D )
  {
    eventPoints.push_back( { p.x(), false } );
    eventPoints.push_back( { p.y(), true  } );
  }

  std::sort( eventPoints.begin(), eventPoints.end() );

  // FIXME: Remove after debugging
  for( auto&& ep : eventPoints )
    std::cout << ep.value << "," << ep.destroyer << "\n"; 
}

} // namespace aleph

#endif
