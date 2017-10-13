#ifndef ALEPH_PERSISTENCE_DIAGRAMS_ENVELOPE_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_ENVELOPE_HH__

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <algorithm>
#include <iterator>
#include <set>
#include <vector>

#include <cassert>

// FIXME: remove after debugging
#include <iostream>

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
  template <class T> void operator()( const PersistenceDiagram<T>& D )
  {
    using Point = typename PersistenceDiagram<T>::Point;

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

      std::cerr << "X = " << it->x() << ", Y = " << it->y() << "\n";
    }
  }
};

} // namespace aleph

#endif
