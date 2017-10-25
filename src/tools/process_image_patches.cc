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
#include <iterator>
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

template <class InputIterator> DataType contrastNorm( InputIterator begin, InputIterator end )
{
  using T = typename std::iterator_traits<InputIterator>::value_type;
  std::vector<T> data( begin, end );

  unsigned n = static_cast<unsigned>( std::sqrt( data.size() ) );

  aleph::math::KahanSummation<T> difference = T();

  for( unsigned i = 0; i < n; i++ )
  {
    for( unsigned j = 0; j < n; j++ )
    {
      auto index = n * i + j;

      if( j+1 <  n )
        difference += std::pow( data.at(index) - data.at(index+1), T(2) );

      if( j   >= 1 )
        difference += std::pow( data.at(index) - data.at(index-1), T(2) );

      if( i+1 <  n )
        difference += std::pow( data.at(index) - data.at(index+n), T(2) );

      if( i-1 <  n )
        difference += std::pow( data.at(index) - data.at(index-n), T(2) );
    }
  }

  return static_cast<DataType>( difference );
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
  // 3. Normalize by the contrast norm

  using IndexType                = decltype( pointCloud.size() );
  PointCloud processedPointCloud = PointCloud( pointCloud.size(), pointCloud.dimension() );

  std::vector<DataType> contrastNorms;
  contrastNorms.reserve( pointCloud.size() );

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

    auto norm = contrastNorm( p.begin(), p.end() );
    contrastNorms.push_back( norm );

    if( norm > DataType() )
    {
      std::transform( p.begin(), p.end(), p.begin(),
                      [&norm] ( DataType x )
                      {
                        return x / norm;
                      } );
    }

    processedPointCloud.set( i, p.begin(), p.end() );
  }

  // Filter patches based on norm --------------------------------------
  //
  // In the original paper, only the top 20% of the contrast norms are
  // being kept. This tool uses a configurable threshold.

  {
    auto contrastNorms_          = contrastNorms;
    double contrastNormThreshold = 0.20; // TODO: make configurable

    std::sort( contrastNorms_.begin(),
               contrastNorms_.end() );

    auto threshold = contrastNorms_.at( std::size_t( std::ceil( (1.0 - contrastNormThreshold) * double( contrastNorms.size() ) ) ) );

    auto numRemainingPatches
      = std::count_if( contrastNorms.begin(), contrastNorms.end(),
                       [&threshold] ( DataType norm )
                       {
                         return norm >= threshold;
                       } );

    PointCloud filteredPointCloud( numRemainingPatches, pointCloud.dimension() );
    IndexType j = 0;

    for( std::size_t i = 0; i < contrastNorms.size(); i++ )
    {
      if( contrastNorms.at(i) >= threshold )
      {
        auto p = processedPointCloud[i];
        filteredPointCloud.set( j++, p.begin(), p.end() );
      }
    }

    swap( processedPointCloud, filteredPointCloud );
  }
}
