#include <aleph/math/SymmetricMatrix.hh>

#include <aleph/persistenceDiagrams/Distances.hh>
#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/BipartiteAdjacencyMatrix.hh>

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include <getopt.h>

// These declarations should remain global because we have to refer to
// them in utility functions that are living outside of `main()`.
using DataType          = double;
using VertexType        = unsigned short;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

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
  std::transform( K.begin(), K.end(), std::back_inserter( simplices ),
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

SimplicialComplex makeLowerFiltration( const SimplicialComplex& K )
{
  auto L = makeFiltration( K );
  L.sort(
    aleph::topology::filtrations::Data<Simplex, std::greater<DataType> >()
  );

  return L;
}

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

int main( int argc, char** argv )
{
  bool absolute              = false;
  bool cycles                = false;
  bool filtration            = false;
  bool minimum               = false;
  bool normalize             = false;
  bool calculateDiagrams     = false;
  bool calculateTrajectories = false;

  {
    static option commandLineOptions[] =
    {
      { "absolute"            , no_argument, nullptr, 'a' },
      { "cycles"              , no_argument, nullptr, 'c' },
      { "filtration"          , no_argument, nullptr, 'f' },
      { "minimum"             , no_argument, nullptr, 'm' },
      { "normalize"           , no_argument, nullptr, 'n' },
      { "persistence-diagrams", no_argument, nullptr, 'p' },
      { "trajectories"        , no_argument, nullptr, 't' },
      { nullptr               , 0          , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "acfmnpt", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'a':
        absolute = true;
        break;
      case 'c':
        cycles = true;
        break;
      case 'f':
        filtration = true;
        break;
      case 'm':
        minimum = true;
        break;
      case 'n':
        normalize = true;
        break;
      case 'p':
        calculateDiagrams = true;
        break;
      case 't':
        calculateTrajectories = true;
        break;
      default:
        break;
      }
    }
  }

  // 1. Read simplicial complexes --------------------------------------

  std::vector<SimplicialComplex> simplicialComplexes;
  simplicialComplexes.reserve( static_cast<unsigned>( argc - optind - 1 ) );

  std::vector<DataType> minData;
  std::vector<DataType> maxData;

  minData.reserve( simplicialComplexes.size() );
  maxData.reserve( simplicialComplexes.size() );

  {
    aleph::topology::io::BipartiteAdjacencyMatrixReader reader;

    if( absolute )
      reader.setUseAbsoluteValues();

    if( minimum )
      reader.setAssignMinimumVertexWeight();

    for( int i = optind; i < argc; i++ )
    {
      auto filename = std::string( argv[i] );

      std::cerr << "* Processing " << filename << "...";

      SimplicialComplex K;
      reader( filename, K );

      std::cerr << "finished\n";

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

      simplicialComplexes.emplace_back( K );
    }
  }

  // 2. Calculate persistent homology ----------------------------------

  using Matrix = aleph::math::SymmetricMatrix<double>;

  // Stores the distance matrix for the trajectories of persistence
  // diagrams. This will only be used if the client set the correct
  // flag.
  Matrix trajectoryDistances;
  if( calculateTrajectories )
    trajectoryDistances = Matrix( simplicialComplexes.size() );

  // Stores the zeroth persistence diagram for calculating trajectories
  // later on. This may need to be extended in order to handle diagrams
  // with higher-dimensional features.
  std::vector<PersistenceDiagram> trajectoryDiagrams;
  if( calculateTrajectories )
    trajectoryDiagrams.reserve( simplicialComplexes.size() );

  for( std::size_t i = 0; i < simplicialComplexes.size(); i++ )
  {
    // The persistence diagram that will be used in the subsequent
    // analysis. This does not necessarily have to stem from data,
    // but can be calculated from a suitable transformation.
    PersistenceDiagram D;

    bool dualize                    = true;
    bool includeAllUnpairedCreators = cycles;

    auto&& K                        = simplicialComplexes[i];

    if( filtration )
    {
      auto L = makeLowerFiltration( K );
      auto U = makeUpperFiltration( K );

      auto lowerDiagrams = aleph::calculatePersistenceDiagrams(
        L,
        dualize,
        includeAllUnpairedCreators
      );

      auto upperDiagrams = aleph::calculatePersistenceDiagrams(
        U,
        dualize,
        includeAllUnpairedCreators
      );

      D = merge( lowerDiagrams.back(), upperDiagrams.back() );
    }
    else
    {
      auto diagrams = aleph::calculatePersistenceDiagrams(
        K,
        dualize,
        includeAllUnpairedCreators
      );

      D = diagrams.back(); // Use the *last* diagram of the filtration so that
                           // we get features in the highest dimension.
    }

    D.removeDiagonal();
    if( !cycles )
      D.removeUnpaired();

    if( normalize )
    {
      // Ensures that points are in [-1:+1] or [0:1] if absolute values
      // have been selected by the client.
      std::transform( D.begin(), D.end(), D.begin(),
        [&absolute, &i, &minData, &maxData] ( const Point& p )
        {
          auto x = p.x();
          auto y = p.y();

          if( absolute )
          {
            x = (x - minData[i]) / (maxData[i] - minData[i]);
            y = (y - minData[i]) / (maxData[i] - minData[i]);
          }
          else
          {
            x = 2 * (x - minData[i]) / (maxData[i] - minData[i]) - 1;
            y = 2 * (y - minData[i]) / (maxData[i] - minData[i]) - 1;
          }

          return Point( x,y );
        }
      );
    }

    if( cycles )
    {
      std::transform( D.begin(), D.end(), D.begin(),
        [] ( const Point& p )
        {
          auto x = p.x();
          auto y = DataType();

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
    // - Calculate persistence diagram trajectories
    // - Calculate 2-norm of the persistence diagrams

    if( calculateDiagrams )
      std::cout << D << "\n\n";
    else if( calculateTrajectories )
      trajectoryDiagrams.push_back( D );
    else
      std::cout << i << "\t" << aleph::pNorm( D ) << "\n";
  }

  // Need to calculate the trajectories afterwards because they require
  // building a database of persistence diagrams.
  if( calculateTrajectories )
  {
    for( std::size_t i = 0; i < trajectoryDiagrams.size(); i++ )
    {
      auto&& Di = trajectoryDiagrams[i];

      for( std::size_t j = i+1; j < trajectoryDiagrams.size(); j++ )
      {
        auto&& Dj = trajectoryDiagrams[j];
        auto dist = aleph::distances::hausdorffDistance(
          Di, Dj
        );

        trajectoryDistances(i,j) = dist;
      }
    }

    // FIXME: replace with proper layout
    std::cout << trajectoryDistances;
  }
}
