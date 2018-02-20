/*
  This is an example file shipped by 'Aleph - A Library for Exploring
  Persistent Homology'.

  This example demonstrates how to calculate a *witness complex* from an
  unstructured point cloud (using Euclidean distances) and calculate its
  persistent homology.

  Demonstrated classes:

    - aleph::PersistenceDiagram
    - aleph::containers::PointCloud
    - aleph::geometry::FLANN
    - aleph::geometry::BruteForce

  Demonstrated functions:

    - aleph::geometry::buildWitnessComplex
    - aleph::calculatePersistenceDiagrams
    - aleph::utilities::convert

  Original author: Bastian Rieck
*/

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/WitnessComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/utilities/String.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <iostream>
#include <string>
#include <vector>

#include <getopt.h>

using namespace aleph;
using namespace containers;
using namespace geometry;
using namespace distances;

void usage()
{
  std::cerr << "Usage: witness_complex FILE [DIMENSION]\n"
            << "\n"
            << "Calculates the witness complex of an unstructured point cloud,\n"
            << "stored in FILE. Euclidean distances are used for the expansion\n"
            << "process. Other optional parameters can be adjusted in order to\n"
            << "change the complex that is built. An optional second argument,\n"
            << "indicating the DIMENSION, can be used to truncate the complex,\n"
            << "making it easier to handle.\n"
            << "\n";
}

int main( int argc, char** argv )
{
  // We first have to specify the data type to use for the subsequent
  // expansion of the witness complex, as well as a distance functor,
  // as this influences the point cloud data type.
  using DataType   = double;
  using PointCloud = PointCloud<DataType>;
  using Distance   = Euclidean<DataType>;

  static option commandLineOptions[] =
  {
    { "landmarks"     , required_argument, nullptr, 'l' },
    { "nu"            , required_argument, nullptr, 'n' },
    { "radius"        , required_argument, nullptr, 'r' },
    { "random"        , no_argument      , nullptr, 'R' },
    { nullptr         , 0                , nullptr,  0  }
  };

  double  landmarksFraction = 0.10; // By default, use 10% of the data points as landmarks
                                    // for the witness complex.
  unsigned nu               = 2;
  DataType radius           = DataType();
  bool randomLandmarks      = false;

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "l:n:r:R", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'l':
        landmarksFraction = std::stod( optarg );
        break;
      case 'n':
        nu = unsigned( std::stoul( optarg ) );
        break;
      case 'r':
        radius = DataType( std::stod( optarg) );
        break;
      case 'R':
        randomLandmarks = true;
        break;
      }
    }
  }

  if( ( argc - optind ) <= 1 )
  {
    usage();
    return -1;
  }

  std::string input = argv[optind++];

  // This loads the point cloud from an unstructured file. The point
  // cloud loader is smart enough to handle things such as different
  // separators in a file.
  auto pointCloud   = aleph::containers::load<DataType>( input );
  auto dimension    = static_cast<unsigned>( pointCloud.dimension() + 1 );
  auto numLandmarks = static_cast<std::size_t>( pointCloud.size() * landmarksFraction );

  if( ( argc - optind ) >= 2 )
    dimension = static_cast<unsigned>( std::stoul( argv[optind++] ) );

  std::vector<std::size_t> landmarks;

  if( randomLandmarks )
  {
    std::cerr << "* Generating landmarks using random strategy...";

    generateRandomLandmarks( pointCloud.size(),
                             numLandmarks,
                             std::back_inserter( landmarks ) );

    std::cerr << "finished\n";
  }
  else
  {
    std::cerr << "* Generating landmarks using max--min strategy...";

    generateMaxMinLandmarks( pointCloud,
                             numLandmarks,
                             std::back_inserter( landmarks ), Distance() );

    std::cerr << "finished\n";
  }

  std::cerr << "* Calculating witness complex with nu=" << nu << ", R=" << radius << ", and d=" << dimension << "...";

  // Aleph gives you some options for detecting optional features such
  // as the wrapper for neighbourhood calculations. I would recommend
  // using the FLANN wrapper. It uses kd-trees, is reasonably fast in
  // higher dimensions and well tested.
  //
  // If this wrapper is not available, we fall back to the brute force
  // wrapper, which---you guessed it---enumerates all neighbours by,
  // well, brute force.
  //
  // If you want to write a wrapper for another library, take a look at
  // the interface of the FLANN wrapper.
  #ifdef ALEPH_WITH_FLANN
    using NearestNeighbours = FLANN<PointCloud, Distance>;
  #else
    using NearestNeighbours = BruteForce<PointCloud, Distance>;
  #endif

  NearestNeighbours nearestNeighbours( pointCloud );

  auto K
    = buildWitnessComplex<Distance>( pointCloud,
                                     landmarks.begin(), landmarks.end(),
                                     dimension,
                                     nu,
                                     radius );

  std::cerr << "finished\n"
            << "* Obtained simplicial complex with " << K.size() << " simplices\n";

  std::cerr << "* Calculating persistence diagrams...";

  // Finally, this function will calculate all persistence diagrams of
  // the simplicial complex. Again, this is a convenience function which
  // assumes that the complex is already in filtration order.
  auto diagrams
    = aleph::calculatePersistenceDiagrams( K );

  std::cerr << "finished\n"
            << "* Obtained " << diagrams.size() << " persistence diagrams\n";

  for( auto&& D : diagrams )
  {
    // Removes all features of zero persistence. They only clutter up
    // the diagonal.
    D.removeDiagonal();

    // This output contains a sort of header (in gnuplot style) so that
    // it is possible to store multiple persistence diagrams in the same
    // file.
    //
    // Note that for the output, it would also be possible just to loop
    // over the individual points of the persistence diagram.
    std::cout << "# Persistence diagram <" << input << ">\n"
              << "#\n"
              << "# Dimension: " << D.dimension() << "\n"
              << "# Entries  : " << D.size() << "\n"
              << D << "\n\n";
  }
}
