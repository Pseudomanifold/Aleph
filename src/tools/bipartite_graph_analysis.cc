#include <aleph/math/SymmetricMatrix.hh>

#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/BipartiteAdjacencyMatrix.hh>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <getopt.h>

// These declarations should remain global because we have to refer to
// them in utility functions that are living outside of `main()`.
using DataType          = double;
using VertexType        = unsigned short;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

SimplicialComplex makeSemiFiltration( const SimplicialComplex& K, bool upper = false )
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

SimplicialComplex makeLowerFiltration( const SimplicialComplex& K, bool reverse = false )
{
  auto L = makeSemiFiltration( K );

  if( reverse )
  {
    L.sort(
      aleph::topology::filtrations::Data<Simplex, std::less<DataType> >()
    );
  }
  else
  {
    L.sort(
      aleph::topology::filtrations::Data<Simplex, std::greater<DataType> >()
    );
  }

  return L;
}

SimplicialComplex makeUpperFiltration( const SimplicialComplex& K, bool reverse = false )
{
  auto L = makeSemiFiltration( K, true );

  if( reverse )
  {
    L.sort(
      aleph::topology::filtrations::Data<Simplex, std::less<DataType> >()
    );
  }
  else
  {
    L.sort(
      aleph::topology::filtrations::Data<Simplex, std::greater<DataType> >()
    );
  }

  return L;
}

SimplicialComplex makeAbsoluteFiltration( const SimplicialComplex& K, bool reverse = false )
{

  auto L = K;

  if( reverse )
  {
    auto functor = [] ( const Simplex& s, const Simplex& t )
    {
      auto w1 = s.data();
      auto w2 = t.data();

      if( std::abs( w1 ) > std::abs( w2 ) )
        return true;
      else if( std::abs( w1 ) == std::abs( w2 ) )
      {
        // This amounts to saying that w1 is negative and w2 is positive,
        // thereby ensuring that the order is consistent.
        if( w1 < w2 )
          return true;
        else
        {
          if( s.dimension() < t.dimension() )
            return true;

          // Absolute value is equal, signed value is equal, and the
          // dimension is equal. We thus have to fall back to merely
          // using the lexicographical order.
          else
            return s < t;
        }
      }

      return false;
    };

    L.sort( functor );
  }
  else
  {
    auto functor = [] ( const Simplex& s, const Simplex& t )
    {
      auto w1 = s.data();
      auto w2 = t.data();

      if( std::abs( w1 ) < std::abs( w2 ) )
        return true;
      else if( std::abs( w1 ) == std::abs( w2 ) )
      {
        // This amounts to saying that w1 is negative and w2 is positive,
        // thereby ensuring that the order is consistent.
        if( w1 < w2 )
          return true;
        else
        {
          if( s.dimension() < t.dimension() )
            return true;

          // Absolute value is equal, signed value is equal, and the
          // dimension is equal. We thus have to fall back to merely
          // using the lexicographical order.
          else
            return s < t;
        }
      }

      return false;
    };

    L.sort( functor );
  }

  return L;
}

using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using Point              = typename PersistenceDiagram::Point;

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

SimplicialComplex applyFiltration( const SimplicialComplex& K,
                                   const std::string& strategy,
                                   bool reverse = false )
{
  auto L = K;

  if( strategy == "standard" )
  {
    if( reverse )
    {
      L.sort(
        aleph::topology::filtrations::Data<Simplex, std::greater<DataType> >()
      );
    }
    else
    {
      L.sort(
        aleph::topology::filtrations::Data<Simplex, std::less<DataType> >()
      );
    }
  }
  else if( strategy == "absolute" )
  {
    if( reverse )
    {
      auto functor = [] ( const Simplex& s, const Simplex& t )
      {
        auto w1 = s.data();
        auto w2 = t.data();

        if( std::abs( w1 ) > std::abs( w2 ) )
          return true;
        else if( std::abs( w1 ) == std::abs( w2 ) )
        {
          // This amounts to saying that w1 is positive and w2 is
          // negative.
          if( w1 > w2 )
            return true;
          else
          {
            if( s.dimension() < t.dimension() )
              return true;

            // Absolute value is equal, signed value is equal, and the
            // dimension is equal. We thus have to fall back to merely
            // using the lexicographical order.
            else
              return s < t;
          }
        }

        return false;
      };

      L.sort( functor );
    }
    else
    {
      auto functor = [] ( const Simplex& s, const Simplex& t )
      {
        auto w1 = s.data();
        auto w2 = t.data();

        if( std::abs( w1 ) < std::abs( w2 ) )
          return true;
        else if( std::abs( w1 ) == std::abs( w2 ) )
        {
          // This amounts to saying that w1 is negative and w2 is
          // positive.
          if( w1 < w2 )
            return true;
          else
          {
            if( s.dimension() < t.dimension() )
              return true;

            // Absolute value is equal, signed value is equal, and the
            // dimension is equal. We thus have to fall back to merely
            // using the lexicographical order.
            else
              return s < t;
          }
        }

        return false;
      };

      L.sort( functor );
    }
  }

  return L;
}

SimplicialComplex assignVertexWeights( const SimplicialComplex& K,
                                       const std::string& strategy,
                                       bool reverse = false )
{
  DataType minData = std::numeric_limits<DataType>::max();
  DataType maxData = std::numeric_limits<DataType>::lowest();

  for( auto&& s : K )
  {
    if( s.dimension() != 1 )
      continue;

    minData = std::min( minData, s.data() );
    maxData = std::max( maxData, s.data() );
  }

  // Setting up the weights --------------------------------------------
  //
  // This function assumes that the simplicial complex is already in
  // filtration ordering with respect to its weights. Hence, we only
  // have to take the *first* weight that we encounter (when using a
  // global vertex weight assignment) or the *extremal* value, which
  // is either a minimum or a maximum depending on the direction.

  std::unordered_map<VertexType, DataType> weight;

  for( auto&& s : K )
  {
    if( s.dimension() != 1 )
      continue;

    auto u     = s[0];
    auto v     = s[1];
    DataType w = DataType(); // weight to assign; depends on filtration

    // Assign the global minimum or maximum. This is rather wasteful
    // because the values do not change, but at least the code makes
    // it clear that all updates are done in the same place.
    if( strategy == "global" )
      w = reverse ? maxData : minData;
    else if( strategy == "local" )
      w = s.data();
    else
      throw std::runtime_error( "Unknown update strategy '" + strategy + "'" );

    // This only performs the update *once*.
    weight.insert( {u,w} );
    weight.insert( {v,w} );
  }

  // Assign the weights ------------------------------------------------
  //
  // Having set up the map of weights, we now only need to traverse it
  // in order to assign weights afterwards.

  auto L = K;

  for( auto it = L.begin(); it != L.end(); ++it )
  {
    if( it->dimension() == 0 )
    {
      auto s = *it;  // simplex
      auto v = s[0]; // vertex

      s.setData( weight.at(v) );

      auto result = L.replace( it, s );
      if( !result )
        throw std::runtime_error( "Unable to replace simplex in simplicial complex" );
    }
  }

  return L;
}

template <class Reader> std::vector<SimplicialComplex> loadSimplicialComplexes( int argc, char** argv, const std::string& minimum )
{
  Reader reader;

  if( minimum == "local" )
    reader.setAssignMinimumVertexWeight();
  else if( minimum == "local_abs" )
    reader.setAssignMinimumAbsoluteVertexWeight();

  std::vector<SimplicialComplex> simplicialComplexes;
  simplicialComplexes.reserve( static_cast<std::size_t>( argc - optind ) );

  for( int i = optind; i < argc; i++ )
  {
    auto filename = std::string( argv[i] );

    std::cerr << "* Processing " << filename << "...";

    SimplicialComplex K;
    reader( filename, K );

    std::cerr << "finished\n";

    simplicialComplexes.emplace_back( K );
  }

  return simplicialComplexes;
}

int main( int argc, char** argv )
{
  bool bipartite             = false;
  bool normalize             = false;
  bool reverse               = false;
  bool verbose               = false;
  bool calculateDiagrams     = false;

  // The default filtration sorts simplices by their weights. Negative
  // weights are treated as being less relevant than positive ones.
  std::string filtration = "standard";

  // Defines how the minimum value for the vertices is to be set. Valid
  // options include:
  //
  //  - global    (uses the global minimum)
  //  - local     (uses the local minimum over all neighbours)
  std::string minimum = "global";

  {
    static option commandLineOptions[] =
    {
      { "bipartite"           , no_argument,       nullptr, 'b' },
      { "normalize"           , no_argument,       nullptr, 'n' },
      { "persistence-diagrams", no_argument,       nullptr, 'p' },
      { "reverse"             , no_argument,       nullptr, 'r' },
      { "verbose"             , no_argument,       nullptr, 'v' },
      { "filtration"          , required_argument, nullptr, 'f' },
      { "minimum"             , required_argument, nullptr, 'm' },
      { nullptr               , 0                , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "bnprtvf:m:", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'b':
        bipartite = true;
        break;
      case 'f':
        filtration = optarg;
        break;
      case 'm':
        minimum = optarg;
        break;
      case 'n':
        normalize = true;
        break;
      case 'p':
        calculateDiagrams = true;
        break;
      case 'r':
        reverse = true;
        break;
      case 'v':
        verbose = true;
        break;
      default:
        break;
      }
    }

    // Check filtration validity ---------------------------------------

    if(    filtration != "absolute"
        && filtration != "standard" )
    {
      std::cerr << "* Invalid filtration value '" << filtration << "', so falling back to standard one\n";
      filtration = "standard";
    }

    // Check minimum validity ------------------------------------------

    if(    minimum != "global"
        && minimum != "local" )
    {
      std::cerr << "* Invalid minimum value '" << minimum << "', so falling back to global one\n";
      minimum = "global";
    }
  }

  // Be verbose about parameters ---------------------------------------

  if( bipartite )
    std::cerr << "* Mode: reading bipartite adjacency matrices\n";
  else
    std::cerr << "* Mode: reading edge lists\n";

  std::cerr << "* Filtration: " << filtration
            << " (" << ( reverse ? "" : "not " ) << "reversed" << ")\n"
            << "* Vertex weight assignment strategy: " << minimum << "\n";

  if( verbose )
    std::cerr << "* Verbose output\n";

  // 1. Read simplicial complexes --------------------------------------

  std::vector<SimplicialComplex> simplicialComplexes;
  simplicialComplexes.reserve( static_cast<unsigned>( argc - optind - 1 ) );

  std::vector<DataType> minData;
  std::vector<DataType> maxData;

  minData.reserve( simplicialComplexes.size() );
  maxData.reserve( simplicialComplexes.size() );

  if( argc - optind >= 1 )
  {
    if( bipartite )
    {
      using Reader = aleph::topology::io::BipartiteAdjacencyMatrixReader;

      simplicialComplexes
        = loadSimplicialComplexes<Reader>( argc, argv, minimum );
    }
  }
  else
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

    for( unsigned i = 0; i < 1e5; i++ )
    {
      auto K
        = makeRandomStratifiedGraph( {2,3}, // FIXME: {2,3,1} for the complete network
                                     engine,
                                     distribution
      );

      simplicialComplexes.emplace_back( K );
    }
  }

  // Determine minimum and maximum values for each complex -------------

  for( auto&& K : simplicialComplexes )
  {
    DataType minData_ = std::numeric_limits<DataType>::max();
    DataType maxData_ = std::numeric_limits<DataType>::lowest();

    // *Always* determine minimum and maximum weights so that we may
    // report them later on. They are only used for normalization in
    // the persistence diagram calculation step.
    for( auto&& s : K )
    {
      minData_ = std::min( minData_, s.data() );
      maxData_ = std::max( maxData_, s.data() );
    }

    minData.push_back( minData_ );
    maxData.push_back( maxData_ );
  }

  // Establish filtration order ----------------------------------------


  // 2. Calculate persistent homology ----------------------------------

  for( std::size_t i = 0; i < simplicialComplexes.size(); i++ )
  {
    // The persistence diagram that will be used in the subsequent
    // analysis. This does not necessarily have to stem from data,
    // but can be calculated from a suitable transformation.
    PersistenceDiagram D;

    auto&& K = simplicialComplexes[i];

    if( filtration == "absolute" )
    {
      auto L        = makeAbsoluteFiltration( K, reverse );
      auto diagrams = aleph::calculatePersistenceDiagrams( L );
      D             = diagrams.back();

      if( verbose )
      {
        std::cerr << "* Absolute value simplicial complex:\n"
                  << L << "\n";
      }
    }
    else if( filtration == "double" )
    {
      auto L = makeLowerFiltration( K, reverse );
      auto U = makeUpperFiltration( K, reverse );

      auto lowerDiagrams = aleph::calculatePersistenceDiagrams( L );
      auto upperDiagrams = aleph::calculatePersistenceDiagrams( U );

      if( !lowerDiagrams.empty() && !upperDiagrams.empty() )
      {
        D = merge(
          lowerDiagrams.back(),
          upperDiagrams.back()
        );
      }

      if( verbose )
      {
        std::cerr << "* Lower simplicial complex:\n"
                  << L << "\n"
                  << "* Upper simplicial complex:\n"
                  << U << "\n";
      }
    }
    else
    {
      if( reverse )
      {
        K.sort(
          aleph::topology::filtrations::Data<Simplex, std::greater<DataType> >()
        );
      }
      else
      {
        K.sort(
          aleph::topology::filtrations::Data<Simplex, std::less<DataType> >()
        );
      }

      if( verbose )
      {
        std::cerr << "* Default simplicial complex:\n"
                  << K << "\n";
      }

      auto diagrams = aleph::calculatePersistenceDiagrams( K );
      D = diagrams.back(); // Use the *last* diagram of the filtration so that
                           // we get features in the highest dimension.
    }

    D.removeDiagonal();
    D.removeUnpaired();

    if( normalize )
    {
      // Ensures that all weights are in [0:1] for the corresponding
      // diagram. This enables the comparison of time-varying graphs
      // or different instances.
      std::transform( D.begin(), D.end(), D.begin(),
        [&i, &minData, &maxData] ( const Point& p )
        {
          auto x = p.x();
          auto y = p.y();

          if( minData[i] != maxData[i] )
          {
            x = (x - minData[i]) / (maxData[i] - minData[i]);
            y = (y - minData[i]) / (maxData[i] - minData[i]);
          }

          return Point( x,y );
        }
      );
    }

    // Determine mode of operation -------------------------------------
    //
    // Several modes of operation exist for this program. They can be
    // set using the flags specified above. At present, the following
    // operations are possible:
    //
    // - Calculate persistence diagrams
    // - Calculate 2-norm of the persistence diagrams

    if( calculateDiagrams )
      std::cout << D << "\n\n";
    else
      std::cout << i << "\t" << aleph::pNorm( D ) << "\n";
  }
}
