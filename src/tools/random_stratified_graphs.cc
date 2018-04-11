#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <chrono>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include <cmath>

using DataType          = float;
using VertexType        = unsigned short;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

template
<
  class Engine,       // random engine to use for weight generation (e.g. std::default_random_engine)
  class Distribution  // distribution to use for weight generation (e.g. std::uniform_real_distribution)
>
SimplicialComplex makeRandomStratifiedGraph(
  const std::vector<unsigned>& strata,
  Engine& engine,
  Distribution& distribution
)
{

  auto n = strata.size();

  if( n <= 1 )
    throw std::runtime_error( "Invalid number of strata" );

  std::vector<Simplex> simplices;

  // Create vertices ---------------------------------------------------
  //
  // The `strata` vector contains the size of each stratum, so we just
  // have to add the correct number of vertices here.

  VertexType index = VertexType(0);
  for( auto&& stratum : strata )
  {
    for( unsigned i = 0; i < stratum; i++ )
      simplices.push_back( Simplex( index++ ) );
  }

  // Create edges ------------------------------------------------------
  //
  // Every stratum is connected to every other stratum, but there are no
  // connections *within* a given stratum.

  VertexType offset = VertexType(0);
  for( decltype(n) i = 0; i < n - 1; i++ )
  {
    // All vertices in the next stratum start with this offset to their
    // indices. It depends on the sum of all vertices in *all* previous
    // strata.
    offset += strata[i];

    for( unsigned j = 0; j < strata[i]; j++ )
    {
      for( unsigned k = 0; k < strata[i+1]; k++ )
      {
        simplices.push_back(
          Simplex(
            {
              VertexType( offset - strata[i] + j ),
              VertexType( offset             + k )
            },
            distribution( engine )
          )
        );
      }
    }
  }

  return SimplicialComplex( simplices.begin(), simplices.end() );
}

int main( int, char** )
{
  std::default_random_engine engine;
  engine.seed(
    static_cast<unsigned>(
      std::chrono::system_clock::now().time_since_epoch().count()
    )
  );

  DataType minWeight = DataType(-1);
  DataType maxWeight = DataType( 1);

  std::uniform_real_distribution<DataType> distribution(
    minWeight,
    std::nextafter( maxWeight, std::numeric_limits<DataType>::max() )
  );

  auto K
    = makeRandomStratifiedGraph( {2,3,1},
                                 engine,
                                 distribution
  );

  std::cerr << K << "\n";
}
