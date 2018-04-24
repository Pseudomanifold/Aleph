/**
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to calculate zero-dimensional persistence diagrams of
  spectra. This is supposed to yield a simple feature descriptor which
  in turn might be used in machine learning pipepline.

  Input:  filename
  Output: persistence diagram

  The persistence diagram represents the superlevel set filtration of
  the input data. This permits us to quantify the number of maxima in
  a data set.

  Original author: Bastian Rieck
*/

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/ConnectedComponents.hh>

#include <aleph/topology/io/FlexSpectrum.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>

#include <cassert>
#include <cmath>

#include <getopt.h>

using DataType           = double;
using VertexType         = unsigned;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

int main( int argc, char** argv )
{
  bool normalize   = false;
  std::string mode = "diagram";

  {
    static option commandLineOptions[] =
    {
      { "mode"     , required_argument, nullptr, 'm' },
      { "normalize", no_argument      , nullptr, 'n' },
      { nullptr    , 0                , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "m:n", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'm':
        mode = optarg;
        break;
      case 'n':
        normalize = true;
        break;
      }
    }
  }

  if( ( argc - optind ) < 1 )
    return -1;

  std::string input = argv[optind++];

  // Parse input -------------------------------------------------------

  std::cerr << "* Reading '" << input << "'...";

  SimplicialComplex K;

  aleph::topology::io::FlexSpectrumReader reader;
  if( normalize )
    reader.normalize( true );

  reader( input, K );

  std::cerr << "finished\n";

  // Calculate persistent homology -------------------------------------

  std::cerr << "* Calculating persistent homology...";

  using PersistencePairing = aleph::PersistencePairing<VertexType>;
  using Traits             = aleph::traits::PersistencePairingCalculation<PersistencePairing>;

  auto&& tuple
    = aleph::calculateZeroDimensionalPersistenceDiagram<Simplex, Traits>( K );

  std::cerr << "finished\n";

  auto&& D       = std::get<0>( tuple );
  auto&& pairing = std::get<1>( tuple );


  // Output ------------------------------------------------------------

  if( mode == "diagram" )
  {
    assert( D.dimension() == 0 );
    assert( D.betti()     == 1 );

    D.removeDiagonal();

    // This ensures that the global maximum is paired with the global
    // minimum of the persistence diagram. This is valid because each
    // function has finite support and is bounded from below.
    std::transform( D.begin(), D.end(), D.begin(),
      [] ( const PersistenceDiagram::Point& p )
      {
        // TODO: we should check whether zero is really the smallest
        // value
        if( p.isUnpaired() )
          return PersistenceDiagram::Point( p.x(), DataType() );
        else
          return PersistenceDiagram::Point( p );
      }
    );

    std::cout << std::setprecision(11) << D << "\n";
  }

  // Transform the (normalized) spectrum into a plane where the
  // $y$-value indicates the persistence of a peak. Thus, it is
  // easier to filter away peaks.
  else if( mode == "transformation" )
  {
    auto&& map = reader.getIndexToValueMap();

    // TODO: there is surely a better way to report the transformed
    // function, but the map sorts values automatically, which sure
    // is nice.
    std::map<double, double> transformedFunction;

    for( auto&& pair : pairing )
    {
      auto&& creator   = pair.first;
      auto&& destroyer = pair.second;
      auto&& sigma     = K.at( creator );
      auto&& tau       = destroyer < K.size() ? K.at( destroyer ) : Simplex( {0,1}, DataType() );

      assert( sigma.dimension() == 0 );
      assert(   tau.dimension() == 1 );

      auto persistence = std::abs( double( sigma.data() ) - double( tau.data() ) );
      auto x           = map.at( sigma[0] );

      transformedFunction[x] = persistence;
    }

    for( auto&& pair : transformedFunction )
      std::cout << std::setprecision(11) << pair.first << "\t" << pair.second << "\n";
  }
}
