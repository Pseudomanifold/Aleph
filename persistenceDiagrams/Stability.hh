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

  bool operator==( const EventPoint& other ) const noexcept
  {
    return value == other.value && destroyer == other.destroyer;
  }

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

  int numActiveFeatures = 0;

  auto numDuplicateValues = [&eventPoints] ( std::size_t i )
  {
    auto eventPoint         = eventPoints.at(i);
    unsigned numOccurrences = 0;
    do
    {
      ++numOccurrences;
      ++i;
    }
    while( i < eventPoints.size() && eventPoints.at(i) == eventPoint );

    return numOccurrences;
  };

  std::vector< std::pair<DataType, DataType> > points;

  for( std::size_t i = 0; i < eventPoints.size(); )
  {
    auto offset = numDuplicateValues(i);

    if( eventPoints.at(i).destroyer )
      numActiveFeatures -= offset;
    else
      numActiveFeatures += offset;

    std::cout << eventPoints.at(i).value << ": " << numActiveFeatures << "\n";

    // FIXME: This is not sanitized; the resulting function will contain
    // duplicate points...
    points.push_back(
      std::make_pair(
        eventPoints.at(i).value,
        numActiveFeatures )
    );

    i += offset;
  }

  return points;
}

} // namespace aleph

#endif
