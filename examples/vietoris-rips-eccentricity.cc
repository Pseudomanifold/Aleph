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

  Demonstrated functions:

    - aleph::geometry::buildVietorisRipsComplex
    - aleph::calculatePersistenceDiagrams
    - aleph::utilities::convert

  Original author: Bastian Rieck

  TODO: Update
*/

#include <aleph/config/FLANN.hh>

#include <aleph/containers/DataDescriptors.hh>
#include <aleph/containers/PointCloud.hh>

#include <aleph/distances/Euclidean.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/utilities/String.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <algorithm>
#include <iostream>
#include <string>

using namespace aleph;
using namespace containers;
using namespace distances;
using namespace geometry;
using namespace utilities;

int main( int argc, char** argv )
{
  if( argc <= 2 )
    return -1;

  using DataType   = double;
  using PointCloud = PointCloud<DataType>;
  using Distance   = Euclidean<DataType>;

  std::string input = argv[1];
  auto pointCloud   = load<DataType>( input );
  auto dimension    = pointCloud.dimension() + 1;
  auto epsilon      = convert<DataType>( argv[2] );
  auto order        = 1u;

  if( argc >= 4 )
    dimension = std::stoul( argv[3] );

  if( argc >= 5 )
    order = static_cast<unsigned>( std::stoul( argv[4] ) );

  // Data descriptor ---------------------------------------------------

  std::cerr << "* Calculating eccentricity data descriptor of order " << order << "...";

  auto eccentricity
    = eccentricities<Distance>( pointCloud, order );

  {
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

#ifdef ALEPH_WITH_FLANN

  FLANN<PointCloud, Distance> flannWrapper( pointCloud );

  auto K
    = buildVietorisRipsComplex( flannWrapper,
                                epsilon,
                                unsigned( dimension ),
                                eccentricity.begin(), eccentricity.end() );

#else

  BruteForce<PointCloud, Distance> flannWrapper( pointCloud );

  auto K
    = buildVietorisRipsComplex( flannWrapper,
                                epsilon,
                                unsigned( dimension ),
                                eccentricity.begin(), eccentricity.end() );

#endif

  std::cerr << "finished\n"
            << "* Obtained simplicial complex with " << K.size() << " simplices\n";

  // Persistent homology -----------------------------------------------

  std::cerr << "* Calculating persistence diagrams...";

  auto diagrams
    = calculatePersistenceDiagrams( K );

  std::cerr << "finished\n"
            << "* Obtained " << diagrams.size() << " persistence diagrams\n";

  for( auto&& D : diagrams )
  {
    D.removeDiagonal();

    std::cout << "# Persistence diagram <" << input << ">\n"
              << "#\n"
              << "# Dimension: " << D.dimension() << "\n"
              << "# Entries  : " << D.size() << "\n"
              << D << "\n\n";
  }
}
