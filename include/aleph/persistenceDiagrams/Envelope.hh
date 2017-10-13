#ifndef ALEPH_PERSISTENCE_DIAGRAMS_ENVELOPE_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_ENVELOPE_HH__

#include <aleph/math/PiecewiseLinearFunction.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <algorithm>
#include <iterator>
#include <set>
#include <vector>

#include <cassert>

namespace aleph
{

/**
  @class Envelope
  @brief Functor for calculating the envelope of a persistence diagram

  The basic idea of this functor is to represent a persistence diagram
  by a simple *envelope function*, i.e. a function that follows maxima
  along the persistence diagram. This can be seen as an easier variant
  of the persistence landscape that does *not* require the calculation
  of intersection points.
*/

class Envelope
{
public:

  /**
    Calculates the *envelope* of a given persistence diagram. i.e. the
    curve connecting the extremal points in the diagram. One can alter
    the functor settings in order to change the way *unpaired* points,
    for example, are being handled.

    @param D Input persistence diagram
    @returns Piecewise linear envelope function
  */

  template <class T> aleph::math::PiecewiseLinearFunction<T> operator()( PersistenceDiagram<T> D )
  {
    using Point = typename PersistenceDiagram<T>::Point;

    if( _removeUnpairedPoints )
      D.removeUnpaired();

    std::vector<Point> P;
    P.reserve( D.size() );

    std::transform( D.begin(), D.end(), std::back_inserter( P ),
                    [] ( const Point& p )
                    {
                      return Point( p.x() + p.y(), p.y() - p.x() );
                    } );

    std::sort( P.begin(), P.end(),
               [] ( const Point& p, const Point& q )
               {
                 if( p.x() != q.x() )
                   return p.x() < q.x();
                 else
                   return p.y() < q.y();
               } );

    std::set<T> X;

    std::transform( P.begin(), P.end(), std::inserter( X, X.begin() ),
                    [] ( const Point& p )
                    {
                      return p.x();
                    } );

    using Coordinate              = std::pair<T, T>;
    using PiecewiseLinearFunction = aleph::math::PiecewiseLinearFunction<T>;

    std::vector<Coordinate> coordinates;

    auto it = P.begin();
    for( auto&& x : X )
    {
      auto next = std::next( it );
      while( next->x() == x )
      {
        ++next;
        ++it;
      }

      assert( x == it->x() );

      coordinates.push_back( std::make_pair( it->x(), it->y() ) );
    }

    return PiecewiseLinearFunction( coordinates.begin(), coordinates.end() );
  }

  void setRemoveUnpairedPoints( bool value = true ) { _removeUnpairedPoints = value; }
  bool    removeUnpairedPoints() const noexcept     { return _removeUnpairedPoints;  }

private:

  /** Flag indicating whether unpaired points should be removed */
  bool _removeUnpairedPoints = true;
};

} // namespace aleph

#endif
