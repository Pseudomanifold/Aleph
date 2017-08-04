#ifndef ALEPH_CONTAINERS_DIMENSIONALITY_ESTIMATORS_HH__
#define ALEPH_CONTAINERS_DIMENSIONALITY_ESTIMATORS_HH__

#include <stdexcept>
#include <vector>

#include <cmath>

#include <aleph/math/KahanSummation.hh>

namespace aleph
{

namespace containers
{

/**
  Estimates local intrinsic dimensionality of a container using its
  nearest neighbours. The underlying assumption of the estimator is
  that points are locally uniformly distributed uniformly. Use this
  estimator with care when analysing unknown data.

  @param container Container to use for dimensionality estimation
  @param k         Number of nearest neighbours
  @param distance  Distance measure

  @returns Vector of local intrinsic dimensionality estimates. Note
           that the numbers are reported *without* rounding.
*/

template <
  class Distance,
  class Container,
  class Wrapper
> std::vector<double> estimateLocalDimensionalityNearestNeighbours( const Container& container,
                                                                    unsigned k,
                                                                    Distance /* distance */ = Distance() )
{
  using IndexType = typename Wrapper::IndexType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<double> > distances;

  Wrapper nnWrapper( container );
  nnWrapper.neighbourSearch( k+1, indices, distances );

  auto n = container.size();

  std::vector<double> estimates;
  estimates.reserve( n );

  for( decltype(n) i = 0; i < n; i++ )
  {
    auto&& nnDistances = distances.at(i);
    auto r1            = aleph::math::accumulate_kahan( nnDistances.begin(), nnDistances.begin() + k,     0.0 ) / static_cast<double>(k  );
    auto r2            = aleph::math::accumulate_kahan( nnDistances.begin(), nnDistances.begin() + k + 1, 0.0 ) / static_cast<double>(k+1);

    estimates.push_back( r1 / ( (r2-r1)*k ) );
  }

  return estimates;
}

/**
  Estimates local intrinsic dimensionality of a container using its
  nearest neighbours. No assumptions about the distribution of data
  points are made. The function uses an iteration over a *range* of
  nearest neighbours and solves a regression problem.

  Please see the publication

    > An evaluation of intrinsic dimensionality estimators\n
    > Peter J. Verveer and Robert P. W. Duin\n
    > IEEE Transactions on Pattern Analysis and Machine Intelligence 17.1, pp. 81-86, 1985

  for more details.

  @param container Container to use for dimensionality estimation

  @param kMin      Minimum number of nearest neighbours to use in
                   computing local dimensionality estimates.

  @param kMax      Maximum number of nearest neighbours to use in
                   computing local dimensionality estimates. This
                   parameter influences performance.

  @param distance  Distance measure

  @returns Vector of local intrinsic dimensionality estimates. Note that
           the numbers are reported *without* rounding.
*/

template <
  class Distance,
  class Container,
  class Wrapper
> std::vector<double> estimateLocalDimensionalityNearestNeighbours( const Container& container,
                                                                    unsigned kMin,
                                                                    unsigned kMax,
                                                                    Distance /* distance */ = Distance() )

{
  if( kMin > kMax )
    std::swap( kMin, kMax );

  if( kMax == 0 || kMin == 0 )
    throw std::runtime_error( "Expecting non-zero number of nearest neighbours" );

  using IndexType   = typename Wrapper::IndexType;
  using ElementType = typename Wrapper::ElementType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<ElementType> > distances;

  Wrapper nnWrapper( container );
  nnWrapper.neighbourSearch( kMax, indices, distances );

  auto n = container.size();

  std::vector<double> estimates;
  estimates.reserve( n );

  for( decltype(n) i = 0; i < n; i++ )
  {
    auto&& nnDistances = distances.at(i);

    std::vector<double> localEstimates;
    localEstimates.reserve( kMax );

    for( unsigned k = kMin; k < kMax; k++ )
    {
      auto r = aleph::math::accumulate_kahan( nnDistances.begin(), nnDistances.begin() + k, 0.0 ) / static_cast<double>(k);
      localEstimates.emplace_back( r );
    }

    // The dimensionality estimates consist of two terms. The first term
    // is similar to the local biased dimensionality estimate.

    std::vector<double> firstTerms;
    firstTerms.reserve( kMax );

    std::vector<double> secondTerms;
    secondTerms.reserve( kMax );

    for( unsigned k = kMin; k < kMax - 1; k++ )
    {
      auto index = k - kMin;
      auto r1    = localEstimates.at(index);
      auto r2    = localEstimates.at(index+1);

      firstTerms.emplace_back ( ( (r2-r1) * r1 ) / k );
      secondTerms.emplace_back( ( (r2-r1) * (r2-r1) ) );
    }

    auto s = aleph::math::accumulate_kahan( firstTerms.begin() , firstTerms.end() , 0.0 );
    auto t = aleph::math::accumulate_kahan( secondTerms.begin(), secondTerms.end(), 0.0 );

    estimates.push_back( s / t );
  }

  return estimates;

}

/**
  Estimates local intrinsic dimensionality of a container using its
  nearest neighbours. No assumptions about the distribution of data
  points are made. The function uses *maximum likelihood estimates*
  for the dimensionality estimates.

  Please see the publication

    > Maximum Likelihood Estimation of Intrinsic Dimension\n
    > Elizaveta Levina and Peter J. Bickel\n
    > Advances in Neural Information Processing Systems, 2005

  for more details.

  @param container Container to use for dimensionality estimation

  @param kMin      Minimum number of nearest neighbours to use in
                   computing local dimensionality estimates.

  @param kMax      Maximum number of nearest neighbours to use in
                   computing local dimensionality estimates. This
                   parameter influences performance.

  @param distance  Distance measure

  @returns Vector of local intrinsic dimensionality estimates. Note that
           the numbers are reported *without* rounding.
*/

template <
  class Distance,
  class Container,
  class Wrapper
> std::vector<double> estimateLocalDimensionalityNearestNeighboursMLE( const Container& container,
                                                                       unsigned kMin,
                                                                       unsigned kMax,
                                                                       Distance /* distance */ = Distance() )
{
  if( kMin > kMax )
    std::swap( kMin, kMax );

  if( kMax == 0 || kMin == 0 )
    throw std::runtime_error( "Expecting non-zero number of nearest neighbours" );

  using IndexType = typename Wrapper::IndexType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<double> > distances;

  Wrapper nnWrapper( container );
  nnWrapper.neighbourSearch( kMax, indices, distances );

  auto n = container.size();

  std::vector<double> estimates;
  estimates.reserve( n );

  for( decltype(n) i = 0; i < n; i++ )
  {
    auto&& nnDistances = distances.at(i);

    std::vector<double> localEstimates;
    localEstimates.reserve( kMax );

    // This follows the notation in the original paper. I dislike using
    // $T_k$ to denote distances, though.
    for( unsigned k = kMin - 1; k < kMax; k++ )
    {
      // Nothing to do here...
      if( k == 0 )
        continue;

      std::vector<double> logEstimates;
      logEstimates.reserve( k-1 );

      for( auto it = nnDistances.begin(); it != nnDistances.begin() + k; ++it )
      {
        if( *it > 0.0 && nnDistances.at(k) > 0.0 )
          logEstimates.push_back( std::log( nnDistances.at(k) / *it ) );

        // This defines log(0) = 0, as usually done in information
        // theory. The original paper does not handle this.
        else
          logEstimates.push_back( 0.0 );
      }

      auto mk = k > 1 ? 1.0 / (k-1) * aleph::math::accumulate_kahan( logEstimates.begin(), logEstimates.end(), 0.0 )
                      : 0.0;

      if( mk > 0.0 )
        mk = 1.0 / mk;
      else
        mk = 0.0;

      localEstimates.push_back( mk );
    }

    estimates.push_back( aleph::math::accumulate_kahan( localEstimates.begin(), localEstimates.end(), 0.0 ) / (kMax - kMin + 1) );
  }

  return estimates;
}

} // namespace containers

} // namespace aleph

#endif
