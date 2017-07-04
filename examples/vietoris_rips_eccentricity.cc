/*
  This is an example file shipped by 'Aleph - A Library for Exploring
  Persistent Homology'.

  This example demonstrates how to obtain a Vietoris--Rips complex
  from an unstructured point cloud (using Euclidean distances) and
  calculate its persistent homology. However, instead of using the
  distances as a source of weights, this program calculates a data
  descriptor based on eccentricities.

  For more information, please refer to the following paper:

    Comparing Dimensionality Reduction Methods Using Data Descriptor Landscapes
    Bastian Rieck and Heike Leitte
    Symposium on Visualization in Data Science (VDS) at IEEE VIS 2015.

  Demonstrated classes:

    - aleph::PersistenceDiagram
    - aleph::containers::PointCloud
    - aleph::geometry::FLANN
    - aleph::geometry::BruteForce
    - aleph::topology::filtrations::LowerStar
    - aleph::topology::filtrations::UpperStar

  Demonstrated functions:

    - aleph::eccentricities
    - aleph::geometry::buildVietorisRipsComplex
    - aleph::calculatePersistenceDiagrams
    - aleph::utilities::convert

  Original author: Bastian Rieck
*/

#include <aleph/config/FLANN.hh>

#include <aleph/containers/DataDescriptors.hh>
#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/utilities/String.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/filtrations/LowerStar.hh>
#include <aleph/topology/filtrations/UpperStar.hh>

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>

void usage()
{
  std::cerr << "Usage: vietoris_rips_eccentricity FILE EPSILON [DIMENSION] [ORDER] [U|L]\n"
            << "\n"
            << "Calculates the Vietoris--Rips complex of an unstructured point\n"
            << "cloud, stored in FILE. Euclidean distances are used during the\n"
            << "expansion process. The maximum distance threshold is specified\n"
            << "by EPSILON. If present, an optional parameter DIMENSION may be\n"
            << "used to truncate the simplicial complex.\n"
            << "\n"
            << "Weights in the simplicial complex will be calculated using the\n"
            << "eccentricity data descriptor. An (optional) ORDER parameter is\n"
            << "used to control how eccentricities are calculated.\n"
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

  auto order = 1u;
  if( argc >= 5 )
    order = static_cast<unsigned>( std::stoul( argv[4] ) );

  // 0 == do not use
  // 1 == upper star
  // 2 == lower star
  unsigned short starFiltration  = 0;
  if( argc >= 6 )
  {
    std::string argument = argv[5];
    argument             = aleph::utilities::trim( argument );
    if( argument == "u" || argument == "U" )
      starFiltration = 1;
    else if( argument == "l" || argument == "L" )
      starFiltration = 2;
  }

  // Data descriptor ---------------------------------------------------
  //
  // The eccentricity data descriptor requires a point cloud and
  // a distance measure for its calculation.

  std::cerr << "* Calculating eccentricity data descriptor of order " << order << "...";

  auto eccentricity
    = aleph::eccentricities<Distance>( pointCloud, order );

  {
    // We transform the values of the data descriptor so that the
    // maximum is mapped to $0$ and the minimum is mapped to $1$,
    // which in essence inverts the values.
    //
    // This is based on a suggestion by Carlsson in his paper
    //
    //  Topological pattern recognition for point cloud data
    //  Acta Numerica, Volume 23, pp. 289--368
    //
    // on page 325. In essence, this ensures that we go through the data
    // from its most outlying points to its inner points. The persistent
    // homology of the data set is likely to change the most during this
    // traversal.

    auto minmax = std::minmax_element( eccentricity.begin(), eccentricity.end() );
    auto min    = *minmax.first;
    auto max    = *minmax.second;

    std::transform( eccentricity.begin(), eccentricity.end(), eccentricity.begin(),
                    [&min, &max] ( double v )
                    {
                      return ( max - v ) / ( max - min );
                    } );
  }

  std::cerr << "finished\n";

  // Vietoris--Rips complex --------------------------------------------

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
    // Since we specified the eccentricity values, the complex will use
    // them as weights. The eccentricities need to be given in the same
    // order as the vertices of the simplicial complex.
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
                                                   unsigned( dimension ),
                                                   eccentricity.begin(), eccentricity.end() );

  #else
    aleph::geometry::BruteForce<PointCloud, Distance> flannWrapper( pointCloud );

    auto K
      = aleph::geometry::buildVietorisRipsComplex( flannWrapper,
                                                   epsilon,
                                                   unsigned( dimension ),
                                                   eccentricity.begin(), eccentricity.end() );
  #endif

  std::cerr << "finished\n"
            << "* Obtained simplicial complex with " << K.size() << " simplices\n";

  switch( starFiltration )
  {
    // We pretend that we do not know the simplex type used by the
    // simplicial complex and get it from its value type instead.
    using Simplex = typename decltype(K)::ValueType;

    case 0:
      break;
    // upper star filtration
    case 1:
    {
      std::cerr << "* Establishing upper-star filtration order...";

      // The upper-star filtration is implemented as a sorting
      // predicate. It uses the eccentricity values calculated
      // above and prepares a lookup table for each simplex.
      aleph::topology::filtrations::UpperStar<Simplex> upperStarFiltration( eccentricity.begin(), eccentricity.end() );

      // Note that for reasons of simplicity, the sorting predicate is
      // always copied during the sorting. To prevent copies, we hence
      // use a reference wrapper from the STL.
      K.sort( std::ref( upperStarFiltration ) );

      std::cerr << "finished\n";
      break;
    }

    // lower-star filtration
    case 2:
    {
      std::cerr << "* Establishing lower-star filtration order...";

      aleph::topology::filtrations::LowerStar<Simplex> lowerStarFiltration( eccentricity.begin(), eccentricity.end() );
      K.sort( std::ref( lowerStarFiltration ) );

      std::cerr << "finished\n";
      break;
    }
  }

  // Persistent homology -----------------------------------------------

  // Finally, this function will calculate all persistence diagrams of
  // the simplicial complex. Again, this is a convenience function which
  // assumes that the complex is already in filtration order.
  std::cerr << "* Calculating persistence diagrams...";

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
