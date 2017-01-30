#ifndef ALEPH_PERSISTENCE_DIAGRAMS_PERSISTENCE_INDICATOR_FUNCTION_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_PERSISTENCE_INDICATOR_FUNCTION_HH__

#include "math/StepFunction.hh"

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <algorithm>
#include <utility>
#include <vector>

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


/**
  Calculates the persistence indicator function of a persistence
  diagram. This function counts the number of 'active' intervals
  for every parameter value. It is a stable summary of a diagram
  and may be used to discern more information about the topology
  and its variation over time.
*/

template <class DataType> aleph::math::StepFunction<DataType> persistenceIndicatorFunction( const PersistenceDiagram<DataType>& D )
{
  using namespace detail;
  using namespace math;

  using EP = EventPoint<DataType>;

  std::vector<EP> eventPoints;
  eventPoints.reserve( 2 * D.size() );

  for( auto&& p : D )
  {
    eventPoints.push_back( { p.x(), false } );
    eventPoints.push_back( { p.y(), true  } );
  }

  std::sort( eventPoints.begin(), eventPoints.end() );

  int numActiveFeatures = 0;

  // Lambda expression for counting duplicate event points that appear after a
  // certain index in the vector of event points. Duplicate event points occur
  // if the persistence diagram contains points with equal values.
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

  StepFunction<DataType> f;

  // Previous interval end point. This is required in order to create proper
  // indicator functions later on.
  DataType previous = DataType();

  for( std::size_t i = 0; i < eventPoints.size(); )
  {
    auto offset = numDuplicateValues(i);

    // Create a new interval if the number of active intervals changes
    if( offset != 0 )
    {
      if( i != 0 && previous != eventPoints.at(i).value )
        f.add( previous, eventPoints.at(i).value, static_cast<DataType>( numActiveFeatures ) );
    }

    if( eventPoints.at(i).destroyer )
      numActiveFeatures -= offset;
    else
      numActiveFeatures += offset;

    // FIXME: This is not sanitized; the resulting function will contain
    // duplicate points...
    points.push_back(
      std::make_pair(
        eventPoints.at(i).value,
        numActiveFeatures )
    );

    previous  = eventPoints.at(i).value;
    i        += offset;
  }

  return f;
}

} // namespace aleph

#endif
