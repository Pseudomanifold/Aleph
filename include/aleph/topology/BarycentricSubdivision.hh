#ifndef ALEPH_TOPOLOGY_BARYCENTRIC_SUBDIVISION_HH__
#define ALEPH_TOPOLOGY_BARYCENTRIC_SUBDIVISION_HH__

#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace aleph
{

namespace topology
{

/**
  A functor for calculating the barycentric subdivision of
  a combinatorial simplicial complex. The interface of the
  functor is just like a regular function call in order to
  permit chaining multiple subdivisions.

  The following code is thus valid:

  \code{.cpp}
  aleph::topology::BarycentricSubdivision f;
  aleph::topology::SimplicialComplex K;

  auto L = f(K);      // first barycentric subdivision
  auto M = f( f(K) ); // second barycentric subdivision

  \endcode
*/

class BarycentricSubdivision
{
public:
  template <class SimplicialComplex> SimplicialComplex operator()( const SimplicialComplex& K ) const
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using VertexType = typename Simplex::VertexType;

    // Stores the new vertex index of the next barycentre vertex that
    // is to be added to the subdivided complex. Initially, this uses
    // the largest vertex index stored in the input data.
    VertexType barycentreVertex = VertexType();

    {
      std::set<VertexType> vertices;
      K.vertices( std::inserter( vertices, vertices.begin() ) );

      if( vertices.empty() )
        return {};

      barycentreVertex = static_cast<VertexType>( *vertices.rbegin() + 1 );
    }

    // Stores the barycentric subdivision of a given simplex. This map
    // contains copies of the 0-simplices. For the other simplices, we
    // store the resulting subdivision as an (unordered) vector.
    std::unordered_map<Simplex, std::vector<Simplex> > subdivision;

    // Stores the subdivided simplicial complex. Initially, this complex
    // only contains the new barycentre vertices. Later on, it will also
    // contain the subdivided skeletons.
    SimplicialComplex L;

    for( auto itDimension = K.begin_dimension(); itDimension != K.end_dimension(); ++itDimension )
    {
      auto&& s = *itDimension;

      if( s.dimension() == 0 )
        subdivision[s] = { s };
      else
      {
        // Copy the data of the old simplex for its barycentric
        // subdivision. Since the subdivision is a _refinement_
        // of the original complex, this makes sense.
        L.push_back( Simplex( barycentreVertex, s.data() ) );

        // Contains all subdivided simplices of the boundary of
        // the current simplex.
        std::vector<Simplex> subdividedBoundary;

        for( auto itBoundary = s.begin_boundary(); itBoundary != s.end_boundary(); ++itBoundary )
        {
          auto pos = K.find( *itBoundary );

          if( pos == K.end() )
            throw std::runtime_error( "Unable to find boundary simplex" );

          auto&& t = *pos;

          subdividedBoundary.insert( subdividedBoundary.end(),
                                     subdivision.at(t).begin(), subdivision.at(t).end() );
        }

        // Cone over the new barycentre vertex and the subdivided
        // boundary of the current simplex. This cone will become
        // the new subdivision of the current simplex.
        std::vector<Simplex> cone;

        for( auto&& t : subdividedBoundary )
        {
          std::vector<VertexType> vertices( t.begin(), t.end() );
          vertices.push_back( barycentreVertex );

          // The new cone simplex needs to use the same weight as the
          // barycentre simplex.
          cone.emplace_back( Simplex( vertices.begin(), vertices.end(),
                                      s.data() ) );
        }

        subdivision[s] = cone;

        // Choose new barycentre vertex for the next simplex. All
        // barycentres will thus be numbered sequentially.
        ++barycentreVertex;
      }
    }

    for( auto&& pair : subdivision )
      L.insert( pair.second.begin(), pair.second.end() );

    return L;
  }
private:

  /**
    Collects all boundaries of the given simplex and returns them. This
    is required in order to update lower-dimensional simplices of cones
    correctly.

    Note that the algorithm is rather stupid and slow; it iterates over
    all faces until the set of all simplices remains fixed.
  */

  template <class Simplex> static std::vector<Simplex> collectBoundaries( const Simplex& s )
  {
    std::set<Simplex> simplices;

    for( auto itBoundary = s.begin_boundary(); itBoundary != s.end_boundary(); ++itBoundary )
      simplices.insert( *itBoundary );

    bool newSimplices = false;
    do
    {
      std::set<Simplex> boundaries;

      for( auto&& s : simplices )
      {
        for( auto itBoundary = s.begin_boundary(); itBoundary != s.end_boundary(); ++itBoundary )
          boundaries.insert( *itBoundary );
      }

      auto oldSize = simplices.size();

      simplices.insert( boundaries.begin(),
                        boundaries.end() );

      auto newSize = simplices.size();
      newSimplices = newSize != oldSize;
    }
    while( newSimplices );


    return { simplices.begin(), simplices.end() };
  }

};

} // namespace topology

} // namespace aleph

#endif
