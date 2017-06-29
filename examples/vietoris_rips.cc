/*
  This is an example file shipped by 'Aleph - A Library for Exploring
  Persistent Homology'.

  This example demonstrates how to obtain a Vietoris--Rips complex
  from an unstructured point cloud (using Euclidean distances) and
  calculate its persistent homology.

  Demonstrated classes:

    - aleph::PersistenceDiagram
    - aleph::containers::PointCloud
    - aleph::geometry::FLANN
    - aleph::geometry::BruteForce

  Demonstrated functions:

    - aleph::geometry::buildVietorisRipsComplex
    - aleph::calculatePersistenceDiagrams
    - aleph::utilities::convert

  Original author: Bastian Rieck
*/

#include <aleph/config/FLANN.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/utilities/String.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <iostream>
#include <string>
#include <vector>

void usage()
{
  std::cerr << "Usage: vietoris_rips FILE EPSILON [DIMENSION]\n"
            << "\n"
            << "Calculates the Vietoris--Rips complex of an unstructured point\n"
            << "cloud, stored in FILE. Euclidean distances are used during the\n"
            << "expansion process. The maximum distance threshold is specified\n"
            << "by EPSILON. If present, an optional parameter DIMENSION may be\n"
            << "used to truncate the simplicial complex.\n"
            << "\n";
}

int main( int argc, char** argv )
{
  if( argc <= 2 )
  {
    usage();
    return -1;
  }

  // We first have to specify the data type to use for the subsequent
  // expansion of the Vietoris--Rips complex. This also results in
  // a different choice of point clouds.
  //
  // Moreover, we specify a distance functor to use for the subsequent
  // expansion process.
  //
  // For educational purposes, we use the full namespace here.
  using DataType   = double;
  using PointCloud = aleph::containers::PointCloud<DataType>;
  using Distance   = aleph::distances::Euclidean<DataType>;

  std::string input = argv[1];

  // This loads the point cloud from an unstructured file. The point
  // cloud loader is smart enough to handle things such as different
  // separators in a file.
  auto pointCloud   = aleph::containers::load<DataType>( input );
  auto dimension    = pointCloud.dimension() + 1;

  // This converts the string supplied by the user to the corresponding
  // data type. Note that the `convert()` function only makes sense for
  // builtin types such as 'double' or 'float'. If you want to use this
  // for your own data types, you need to overload `operator>>` because
  // the converter internally uses `std::stringstream` for tokens.
  auto epsilon = aleph::utilities::convert<DataType>( argv[2] );

  if( argc >= 4 )
    dimension = std::stoul( argv[3] );

  std::cerr << "* Calculating Vietoris--Rips complex with eps=" << epsilon << " and d=" << dimension << "...";

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
    aleph::geometry::FLANN<PointCloud, Distance> flannWrapper( pointCloud );

    // That's really all there is to is: the convenience function below
    // uses a neighbourhood wrapper and additional parameters and
    // creates an appropriate Vietoris--Rips complex.
    //
    // Since we did not specify anything else, the complex will contain
    // the high-dimensional distances as weights. The vertices of K are
    // assigned a value of 0, while the edges are assigned a value
    // according to the distance between their corresponding vertices.
    //
    // The simplicial complex is ordered according to this weight, from
    // low values to high values, resulting in a filtration of the
    // sublevel sets of the distance function.
    //
    // Other filtrations are possible but need to be applied manually
    // afterwards.
    auto K
      = aleph::geometry::buildVietorisRipsComplex( flannWrapper,
                                                   epsilon,
                                                   unsigned( dimension ) );
  #else
    aleph::geometry::BruteForce<PointCloud, Distance> bruteForceWrapper( pointCloud );

    auto K
      = aleph::geometry::buildVietorisRipsComplex( bruteForceWrapper,
                                                   epsilon,
                                                   unsigned( dimension ) );
  #endif

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
