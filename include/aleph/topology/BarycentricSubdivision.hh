#ifndef ALEPH_TOPOLOGY_BARYCENTRIC_SUBDIVISION_HH__
#define ALEPH_TOPOLOGY_BARYCENTRIC_SUBDIVISION_HH__

#include <aleph/utilities/EmptyFunctor.hh>

#include <iterator>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

// Since the data type of the simplex class is allowed to be a boolean
// as well, there will be a multiplication involved in the functor and
// the compiler rightfully warns about this because two boolean values
// are part of the calculation.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"

#ifdef __clang__
  // Nothing to do here, but clang also defines __GNUCC__, making this
  // check below moot.
#else
  #if __GNUC__ < 6
    #pragma GCC diagnostic ignored "-Werror"
  #else
    #pragma GCC diagnostic ignored "-Wint-in-bool-context"
  #endif
#endif

namespace aleph
{

namespace topology
{

/**
  @class BarycentricSubdivision
  @brief Barycentric subdivision functor for complexes

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

  The main operation of the functor is handled by BarycentricSubdivision::operator()(). Please
  refer to the documentation of this function for more information.
*/

class BarycentricSubdivision
{
public:

  /**
    Performs a barycentric subdivision of the given simplicial complex
    and returns the result.

    @param K       Simplicial complex
    @param functor Functor for assigning weights to subdivided
                   simplices. The premise of the functor is to
                   return a factor that determines how the old
                   weight of a given simplex is being used for
                   the subdivided one. The functor hence needs
                   to return a scaling factor. By default, the
                   scaling factor is 1, so the weights will be
                   merely copied.\n

                   The functor needs to provide the evaluation
                   operator `operator()()` taking one unsigned
                   value that describes the *dimension* of the
                   current simplex. It should return a scaling
                   factor for this dimension.\n

                   The subsequent code provides a functor that
                   is a good choice when subdividing edges for
                   which the weight specifies a *length*:\n

                   \code{.cpp}

                    aleph::topology::BarycentricSubdivision subdivision;
                    aleph::topology::SimplicialComplex K;

                    auto L = subdivision( K, [] ( std::size_t dimension ) { return dimension == 0 ? 0 : 0.5; } );
                    L.recalculateWeights( true, true );
                   \endcode\n

                    This will divide the length of every edge,
                    while setting all other weights to zero so
                    that the resulting simplicial complex does
                    model a distance function correctly.

    @tparam SimplicialComplex Simplicial complex class type
    @tparam Functor           Functor type for assigning weights
  */

  template <class SimplicialComplex, class Functor = aleph::utilities::EmptyFunctor> SimplicialComplex operator()( const SimplicialComplex& K, Functor&& functor = Functor() ) const
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using VertexType = typename Simplex::VertexType;
    using DataType   = typename Simplex::DataType;

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
        {
          auto data = s.data() * DataType( functor( 0 ) );
          L.push_back( Simplex( barycentreVertex, data ) );
        }

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
                                      s.data() * DataType( functor( vertices.size() ) ) ) );
        }

       // Cone boundaries; required in order to ensure consistency of
       // the resulting simplicial complex.
        std::vector<Simplex> boundaries
          = collectBoundaries( subdividedBoundary.begin(), subdividedBoundary.end() );

        for( auto&& t : boundaries )
        {
          std::vector<VertexType> vertices( t.begin(), t.end() );
          vertices.push_back( barycentreVertex );

          // The new simplex is directly inserted into the boundary of
          // the resulting simplicial complex---it will not be used by
          // any other simplex during the subdivision.
          L.push_back( Simplex( vertices.begin(), vertices.end(),
                                s.data() * DataType( functor( vertices.size() ) ) ) );
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

  /**
    Overloaded variant of the function above. It permits the collection
    of all boundaries of a range of simplices.
  */

  template <class InputIterator> static auto collectBoundaries( InputIterator begin, InputIterator end ) -> std::vector<typename std::iterator_traits<InputIterator>::value_type>
  {
    using Simplex = typename std::iterator_traits<InputIterator>::value_type;
    std::set<Simplex> simplices;

    for( auto it = begin; it != end; ++it )
    {
      auto&& boundaries = collectBoundaries( *it );

      simplices.insert( boundaries.begin(),
                        boundaries.end() );
    }

    return { simplices.begin(), simplices.end() };
  }
};

} // namespace topology

} // namespace aleph

#pragma GCC diagnostic pop

#endif
