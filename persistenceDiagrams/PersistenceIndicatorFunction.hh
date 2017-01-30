#ifndef ALEPH_PERSISTENCE_DIAGRAMS_PERSISTENCE_INDICATOR_FUNCTION_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_PERSISTENCE_INDICATOR_FUNCTION_HH__

#include "math/StepFunction.hh"

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include <cmath>

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

template <class T> T next( T x )
{
  if( std::numeric_limits<T>::is_integer )
    return x+1;
  else
    return std::nextafter( x, std::numeric_limits<T>::max() );
}

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

  unsigned numActiveFeatures = 0;
  bool isDegenerate          = false;

  // Sanity check: A destroyer and a creator should never have the same
  // function value. Else, the number of 'active' intervals will behave
  // somewhat strangely.
  if( eventPoints.size() >= 2 )
  {
    for( std::size_t i = 0; i < eventPoints.size(); i += 2 )
    {
      if( eventPoints[i].value == eventPoints[i+1].value )
      {
        if( eventPoints[i].destroyer != eventPoints[i+1].destroyer )
        {
          isDegenerate = true;
          break;
        }
      }
    }
  }

  // Lambda expression for counting duplicate event points that appear after a
  // certain index in the vector of event points. Duplicate event points occur
  // if the persistence diagram contains points with equal values.
  auto numDuplicateValues = [&eventPoints] ( std::size_t i )
  {
    auto eventPoint         = eventPoints.at(i);
    unsigned numOccurrences = 0;
    unsigned numDestroyers  = 0;
    do
    {
      numDestroyers += eventPoints.at(i).destroyer;

      ++numOccurrences;
      ++i;
    }
    while( i < eventPoints.size() && eventPoints.at(i).value == eventPoint.value );

    return std::make_pair( numOccurrences, numDestroyers );
  };

  StepFunction<DataType> f;

  // Previous interval end point. This is required in order to create proper
  // indicator functions later on.
  DataType previous = DataType();

  for( std::size_t i = 0; i < eventPoints.size(); )
  {
    auto pair        = numDuplicateValues(i);
    auto occurrences = pair.first;
    auto destroyers  = pair.second;
    auto creators    = occurrences - destroyers;

    bool useNextPoint = false;

    // Case 1: No duplicates or duplicates of the same type. In this case, the
    // number of active intervals changes and we add an interval. It comprises
    // the number of active intervals up to this point.
    if( occurrences == 1 || creators == occurrences || creators == 0 )
    {
      if( i != 0 )
        f.add( previous, eventPoints.at(i).value, static_cast<DataType>( numActiveFeatures ) );
    }
    
    // Case 2: There are duplicate creation & destruction values. This
    // necessitates the creation of two intervals: one interval at the
    // current even point, with the proper number of active intervals,
    // the other one *directly* afterwards to indicate the destruction
    // implied by the values.
    else
    {
      f.add( previous, eventPoints.at(i).value, static_cast<DataType>( numActiveFeatures + creators ) );
      useNextPoint = true;
    }

    numActiveFeatures += creators;
    numActiveFeatures -= destroyers;

    previous  = useNextPoint ? next( eventPoints.at(i).value ) : eventPoints.at(i).value;
    i        += occurrences;
  }

  return f;
}

} // namespace aleph

#endif
