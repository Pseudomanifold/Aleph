#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <algorithm>
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

// FIXME: copied from bipartite graph analysis tool
SimplicialComplex makeFiltration( const SimplicialComplex& K, bool upper = false )
{
  std::vector<Simplex> simplices;
  simplices.reserve( K.size() );

  // Keep vertices if they are above/below the desired data type
  // threshold for the filtration.
  std::transform( K.begin(), K.end(), std::back_inserter( simplices ),
    [&upper] ( const Simplex& s )
    {
      if( ( upper && s.data() > DataType() ) || ( !upper && s.data() < DataType() ) )
        return s;

      // Copy the simplex but set its weight to be zero because it does
      // not correspond to any structure that we want to learn.
      else
      {
        std::vector<VertexType> vertices( s.begin(), s.end() );
        return Simplex( s.begin(), s.end(), DataType() );
      }
    }
  );

  // Ensure that all vertices are created at threshold zero. This
  // indicates that vertices are always available in the network,
  // regardless of weight threshold.
  //
  // FIXME: this somewhat interferes with the weight selection in
  // the reader class; not sure how to merge those aspects
  std::transform( simplices.begin(), simplices.end(), simplices.begin(),
    [] ( const Simplex& s )
    {
      if( s.dimension() == 0 )
        return Simplex( *s.begin(), DataType() );
      else
        return s;
    }
  );

  // Remove higher-dimensional simplices (edges) that do not have
  // a part in the current filtration.
  simplices.erase(
    std::remove_if( simplices.begin(), simplices.end(),
      [] ( const Simplex& s )
      {
        return s.dimension() != 0 && s.data() == DataType();
      }
    ), simplices.end()
  );

  return SimplicialComplex( simplices.begin(), simplices.end() );
}

// FIXME: copied from bipartite graph analysis tool
SimplicialComplex makeLowerFiltration( const SimplicialComplex& K )
{
  auto L = makeFiltration( K );
  L.sort(
    aleph::topology::filtrations::Data<Simplex, std::greater<DataType> >()
  );

  return L;
}

// FIXME: copied from bipartite graph analysis tool
SimplicialComplex makeUpperFiltration( const SimplicialComplex& K )
{
  auto L = makeFiltration( K, true );
  L.sort(
    aleph::topology::filtrations::Data<Simplex, std::less<DataType> >()
  );

  return makeFiltration( K, true );
}

using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using Point              = typename PersistenceDiagram::Point;

// FIXME: copied from bipartite graph analysis tool
PersistenceDiagram merge( const PersistenceDiagram& D, const PersistenceDiagram& E )
{
  PersistenceDiagram F;

  if( D.dimension() != F.dimension() )
    throw std::runtime_error( "Persistence diagram dimensions have to agree" );

  for( auto&& diagram : { D, E } )
    for( auto&& p : diagram )
      F.add( p.x(), p.y() );

  return F;
}

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

  // TODO: make configurable
  bool normalize = true;

  std::uniform_real_distribution<DataType> distribution(
    minWeight,
    std::nextafter( maxWeight, std::numeric_limits<DataType>::max() )
  );

  // TODO: make configurable
  for( unsigned i = 0; i < 1e5; i++ )
  {
    auto K
      = makeRandomStratifiedGraph( {2,3}, // FIXME: {2,3,1} for the complete network
                                   engine,
                                   distribution
    );

    DataType minData = std::numeric_limits<DataType>::max();
    DataType maxData = std::numeric_limits<DataType>::lowest();

    for( auto&& s : K )
    {
      minData = std::min( minData, s.data() );
      maxData = std::max( maxData, s.data() );
    }

    PersistenceDiagram D;

    // This uses the upper--lower filtration, which is not theoretically
    // justified. This should be configurable.
    {
      auto L = makeLowerFiltration( K );
      auto U = makeUpperFiltration( K );

      auto lowerDiagrams = aleph::calculatePersistenceDiagrams( L );
      auto upperDiagrams = aleph::calculatePersistenceDiagrams( U );

      if( !lowerDiagrams.empty() && !upperDiagrams.empty() )
      {
        D = merge(
          lowerDiagrams.front(),
          upperDiagrams.front()
        );
      }
    }

    if( normalize && minData != maxData )
    {
      std::transform( D.begin(), D.end(), D.begin(),
        [&minData, &maxData] ( const Point& p )
        {
          auto x = p.x();
          auto y = p.y();

          x = (x - minData) / (maxData - minData);
          y = (y - minData) / (maxData - minData);

          return Point( x,y );
        }
      );
    }

    D.removeDiagonal();
    D.removeUnpaired();

    if( !D.empty() )
      std::cout << aleph::pNorm( D ) << "\n";
  }
}
