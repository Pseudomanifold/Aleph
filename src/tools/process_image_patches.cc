/**
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Given a set of patches from a database of images, this tool performs
  the pre-processing steps of the corresponding point cloud. Mainly, a
  procedure from the paper *On the Local Behavior of Spaces of Natural
  Images* by Gunnar Carlsson et al. is followed.
*/

#include <aleph/containers/PointCloud.hh>

#include <aleph/math/KahanSummation.hh>

#include <algorithm>
#include <iostream>
#include <string>

#include <cmath>

using DataType   = float;
using PointCloud = aleph::containers::PointCloud<DataType>;

template <class T> T log( T x )
{
  if( x == T() )
    return T();
  else
    return std::log10( x );
}

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  // Input -------------------------------------------------------------
  //
  // This tool assumes that the input is already in the form of a point
  // cloud, containing the 'raw' image patches.

  std::string filename = argv[1];

  std::cerr << "* Loading input point cloud...";

  PointCloud pointCloud = aleph::containers::load<DataType>( filename );

  std::cerr << "finished\n";

  // Pre-processing ----------------------------------------------------
  //
  // 1. Replace values by their logarithm
  // 2. Subtract mean

  using IndexType                = decltype( pointCloud.size() );
  PointCloud processedPointCloud = PointCloud( pointCloud.size(), pointCloud.dimension() );

  for( IndexType i = 0; i < pointCloud.size(); i++ )
  {
    auto p = pointCloud[i];

    std::transform( p.begin(), p.end(), p.begin(),
                    [] ( DataType x )
                    {
                      return log( x );
                    } );

    auto mean  = aleph::math::accumulate_kahan_sorted( p.begin(), p.end(), DataType() );
    mean      /= static_cast<DataType>( pointCloud.dimension() );

    std::transform( p.begin(), p.end(), p.begin(),
                    [&mean] ( DataType x )
                    {
                      return x - mean;
                    } );

    processedPointCloud.set( i, p.begin(), p.end() );
  }

}
