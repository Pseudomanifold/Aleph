#include <aleph/containers/PointCloud.hh>

#include <aleph/config/FLANN.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>

#include <aleph/geometry/distances/Euclidean.hh>
#include <aleph/geometry/distances/Manhattan.hh>
#include <aleph/geometry/distances/Traits.hh>

#include <aleph/math/KahanSummation.hh>

#include <iostream>
#include <string>
#include <vector>

#include <getopt.h>

using DataType          = double;
using PointCloud        = aleph::containers::PointCloud<DataType>;
using EuclideanDistance = aleph::distances::Euclidean<DataType>;
using ManhattanDistance = aleph::distances::Manhattan<DataType>;

template <class Distance> std::vector<DataType> pairwiseDistances( const PointCloud& pointCloud, Distance distance = Distance() )
{
  if( pointCloud.empty() )
    return {};

  std::vector<DataType> distances;
  distances.reserve( pointCloud.size() * ( pointCloud.size() - 1 ) / 2 );

  auto d       = pointCloud.dimension();
  auto n       = pointCloud.size();
  using Traits = aleph::distances::Traits<Distance>;

  Traits traits;
  for( decltype(n) i = 0; i < n; i++ )
  {
    for( decltype(n) j = i+1; j < n; j++ )
    {
      auto dist = distance( pointCloud[i].begin(),
                            pointCloud[j].begin(),
                            d );

      // I want to be sure that I get the *square* of the distance, but
      // I cannot take this transformation within the distance functor,
      // such as the Manhattan distance, for granted.
      dist      = traits.from( dist );
      dist     *= dist;

      distances.emplace_back( dist );
    }
  }

  return distances;
}

template <class Distance> void calculateDegrees( const PointCloud& pointCloud, DataType epsilon, Distance /* distance = Distance() */,
                                                 std::vector<unsigned>& unweightedDegrees,
                                                 std::vector<DataType>& weightedDegrees )
{
  if( epsilon == DataType() )
    return;

#ifdef ALEPH_WITH_FLANN
  using NearestNeighbours = aleph::geometry::FLANN<PointCloud, Distance>;
#else
  using NearestNeighbours = aleph::geometry::BruteForce<PointCloud, Distance>;
#endif

  NearestNeighbours nn( pointCloud );

  using ElementType = typename NearestNeighbours::ElementType;
  using IndexType   = typename NearestNeighbours::IndexType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<ElementType> > distances;

  nn.radiusSearch( epsilon, indices, distances );

  auto n = pointCloud.size();

  unweightedDegrees.clear();
  unweightedDegrees.reserve( n );

  weightedDegrees.clear();
  weightedDegrees.reserve( n );

  for( decltype(n) i = 0; i < n; i++ )
  {
    auto unweighted = static_cast<unsigned>( indices.at(i).size() );
    auto&& D        = distances.at(i);
    auto weighted   = static_cast<DataType>( aleph::math::accumulate_kahan_sorted( D.begin(), D.end(), DataType() ) );

    unweightedDegrees.emplace_back( unweighted );
    weightedDegrees.emplace_back( weighted );
  }
}

template <class Container> void containerAsJSON( std::ostream& out, const Container& container, const std::string& name, unsigned indent = 2 )
{
  if( container.empty() )
    return;

  out << std::string( indent, ' ' ) << "\"" << name << "\": "
      << "[";

  for( auto it = container.begin(); it != container.end(); ++it )
  {
    if( it != container.begin() )
      out << ",";

    out << *it;
  }

  out << "]\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "distance"      , required_argument, nullptr, 'd' },
    { "epsilon"       , required_argument, nullptr, 'e' },
    { nullptr         , 0                , nullptr,  0  }
  };

  std::string selectedDistanceFunctor = "euclidean";
  DataType epsilon                    = DataType();

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "d:e:", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'd':
        selectedDistanceFunctor = optarg;
        break;
      case 'e':
        epsilon = static_cast<DataType>( std::stod( optarg ) );
        break;
      }
    }
  }

  if( (argc - optind) < 1 )
    return -1;

  auto filename   = std::string( argv[optind++] );
  auto pointCloud = aleph::containers::load<DataType>( filename );

  std::cerr << "* Loaded point cloud with " << pointCloud.size() << " points\n";

  std::vector<DataType> distances;
  std::vector<unsigned> unweightedDegrees;
  std::vector<DataType> weightedDegrees;

  if( selectedDistanceFunctor == "euclidean" )
  {
    distances = pairwiseDistances( pointCloud, EuclideanDistance() );

    calculateDegrees( pointCloud, epsilon, EuclideanDistance(),
                      unweightedDegrees, weightedDegrees );
  }
  else if( selectedDistanceFunctor == "manhattan" )
  {
    distances = pairwiseDistances( pointCloud, ManhattanDistance() );

    calculateDegrees( pointCloud, epsilon, ManhattanDistance(),
                      unweightedDegrees, weightedDegrees );
  }

  std::cout << "{\n";

  containerAsJSON( std::cout, distances, "distances" );

  if( !unweightedDegrees.empty() )
  {
    std::cout << ",\n";
    containerAsJSON( std::cout, unweightedDegrees, "unweighted_degrees" );
  }

  if( !weightedDegrees.empty() )
  {
    std::cout << ",\n";
    containerAsJSON( std::cout, weightedDegrees, "weighted_degrees" );
  }

  std::cout << "}\n";
}
