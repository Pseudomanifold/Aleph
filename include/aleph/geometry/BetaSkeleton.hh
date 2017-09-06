#ifndef ALEPH_GEOMETRY_BETA_SKELETON_HH__
#define ALEPH_GEOMETRY_BETA_SKELETON_HH__

#include <aleph/geometry/distances/Traits.hh>

#include <aleph/utilities/ContainerOperators.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <vector>

namespace aleph
{

namespace geometry
{

namespace detail
{

/*
  Describes a d-dimensional ball with a certain diameter and a certain
  centre; the ball will be used during beta-skeleton construction.
*/

class BetaBall
{
public:

  // This is required because the lune cannot initialize its balls prior
  // to calculating certain distances in the point cloud.
  BetaBall()
  {
  }

  template <class T> BetaBall( const std::vector<T>& centre, double diameter )
    : _centre( centre.begin(), centre.end() )
    , _radius( 0.5 * diameter )
  {
  }

  template <class T > bool contains( const std::vector<T>& other ) const
  {
    double distance = 0.0;

    for( std::size_t i = 0; i < _centre.size(); i++ )
      distance += ( _centre[i] - other[i] ) * ( _centre[i] - other[i] );

    return distance <= _radius * _radius;
  }

private:
  std::vector<double> _centre;
  double _radius = 0.0;
};

/*
  Describes the lune that is used as the empty region of the beta-skeleton
  and offers a way of checking for the presence of data points.
*/

template <class Container, class Index = std::size_t> class BetaLune
{
public:
  BetaLune( const Container& container,
            Index p,
            Index q,
            double beta,
            double d )
    : _container( container )
  {
    using namespace aleph::utilities;

    auto&& P      = container[p];
    auto&& Q      = container[q];
    auto centreP  = (1.0-0.5*beta) * P + 0.5*beta * Q;
    auto centreQ  = (1.0-0.5*beta) * Q + 0.5*beta * P;
    auto diameter =  beta * d;

    _pBall = BetaBall( centreP, diameter );
    _qBall = BetaBall( centreQ, diameter );
  }

  bool contains( Index r ) const
  {
    return _pBall.contains( _container[r] )
        && _qBall.contains( _container[r] );
  }

private:
  const Container& _container;

  BetaBall _pBall;
  BetaBall _qBall;
};

} // namespace detail

/**
  Builds a \f$\beta\f$-skeleton for a given container. This skeleton is
  defined as an empty region graph where the empty region is defined by
  two congruent disks, whose diameter is initially the distance between
  two points, scaled by \f$\beta\f$. Edges will be created only if this
  region is devoid of any other points.

  @param container Container from which to calculate the skeleton
  @param beta      Scaling parameter for the empty region
  @param distance  Distance functor to use for the calculation

  @returns Simplicial complex representing the \f$\beta\f$-skeleton.
*/

template <class Distance, class Container, class Index = std::size_t>
  auto buildBetaSkeletonNaive( const Container& container,
                               double beta,
                               Distance distance = Distance() )
    -> topology::SimplicialComplex< topology::Simplex<typename Distance::ResultType, Index> >

{
  using Traits            = aleph::distances::Traits<Distance>;
  using DataType          = typename Distance::ResultType;
  using VertexType        = Index;
  using Simplex           = topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = topology::SimplicialComplex<Simplex>;

  auto n = container.size();
  auto d = container.dimension();

  Traits traits;
  SimplicialComplex betaSkeleton;

  for( Index i = 0; i < n; i++ )
    betaSkeleton.push_back( Simplex(i) );

  for( Index i = 0; i < n; i++ )
  {
    auto&& p = container[i];

    for( Index j = i+1; j < container.size(); j++ )
    {
      auto&& q  = container[j];
      auto dist = traits.from( distance( p.begin(), q.begin(), d ) );

      using namespace detail;

      BetaLune<Container> lune( container,
                                i,
                                j,
                                beta,
                                dist );

      bool addEdge = true;

      for( Index r = 0; r < container.size(); r++ )
      {
        if( r != i && r != j && lune.contains( r ) )
        {
          addEdge = false;
          break;
        }
      }

      if( addEdge )
        betaSkeleton.push_back( Simplex( {i,j}, dist ) );
    }
  }

  return betaSkeleton;
}

} // namespace geometry

} // namespace aleph

#endif
