#include <aleph/topology/io/LexicographicTriangulation.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/PhiPersistence.hh>

#include <aleph/topology/BarycentricSubdivision.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>
#include <aleph/topology/Skeleton.hh>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using DataType          = bool;
using VertexType        = unsigned short;

using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

/**
  Models a signature consisting of Betti numbers, i.e. a set of natural
  numbers. Signatures are comparable and are ordered lexicographically,
  i.e. in the manner one would expect.
*/

class Signature
{
public:
  template <class InputIterator> Signature( InputIterator begin, InputIterator end )
    : _betti( begin, end )
  {
  }

  bool operator<( const Signature& other ) const noexcept
  {
    return std::lexicographical_compare( _betti.begin(), _betti.end(),
                                         other._betti.begin(), other._betti.end() );
  }

  using const_iterator = typename std::vector<std::size_t>::const_iterator;

  const_iterator begin() const noexcept { return _betti.begin(); }
  const_iterator end()   const noexcept { return _betti.end() ;  }

private:
  std::vector<std::size_t> _betti;
};

std::ostream& operator<<( std::ostream& o, const Signature& s )
{
  o << "(";

  for( auto it = s.begin(); it != s.end(); ++it )
  {
    if( it != s.begin() )
      o << ",";

    o << *it;
  }

  o << ")";

  return o;
}

/**
  Converts a vector of persistence diagrams into a Betti signature. The
  maximum dimension needs to be specified in order to ensure that empty
  or missing persistence diagrams can still be handled correctly.
*/

template <class PersistenceDiagram> Signature makeSignature( const std::vector<PersistenceDiagram>& diagrams, std::size_t D )
{
  std::vector<std::size_t> bettiNumbers( D+1 );

  for( auto&& diagram : diagrams )
  {
    auto d = diagram.dimension();
    auto b = diagram.betti();

    bettiNumbers.at(d) = b;
  }

  return Signature( bettiNumbers.begin(), bettiNumbers.end() );
}

/**
  Enumerates all possible perversities for a given dimension. One could
  say that this function is as wicked as possible.
*/

std::vector<aleph::Perversity> getPerversities( unsigned dimension )
{
  std::vector<aleph::Perversity> perversities;

  std::map< unsigned, std::vector<int> > possibleValues;

  for( unsigned d = 0; d < dimension; d++ )
  {
    // Note that no shift in dimensions is required: as the dimension is
    // zero-based, the maximum value of the perversity in dimension zero
    // is zero. This is identical to demanding
    //
    //   -1 \leq p_k \leq k-1
    //
    // for k = [1,...,d].
    for( int k = -1; k <= static_cast<int>( d ); k++ )
      possibleValues[d].push_back( k );
  }

  std::size_t numCombinations = 1;

  for( auto&& pair : possibleValues )
    numCombinations *= pair.second.size();

  // Stores the current index that was reached while traversing all
  // possible values of the perversity. The idea is that to add one
  // to the last index upon adding a new value. If the index is too
  // large, it goes back to zero and the previous index is modified
  // by one.
  std::vector<std::size_t> indices( possibleValues.size() );

  for( std::size_t i = 0; i < numCombinations; i++ )
  {
    std::vector<int> values;

    for( unsigned d = 0; d < dimension; d++ )
    {
      auto index = indices.at(d);
      values.emplace_back( possibleValues.at(d).at(index) );
    }

    // Always increase the last used index by one. If an overflow
    // occurs, the previous indices are updated as well.
    {
      auto last        = dimension - 1;
      indices.at(last) = indices.at(last) + 1;

      // Check & propagate potential overflows to the previous indices
      // in the range. The index vector is reset to zero only when all
      // combinations have been visited.
      if( indices.at(last) >= possibleValues.at(last).size() )
      {
        indices.at(last) = 0;
        for( unsigned d = 1; d < dimension; d++ )
        {
          auto D        = dimension - 1 - d;
          indices.at(D) = indices.at(D) + 1;

          if( indices.at(D) < possibleValues.at(D).size() )
            break;
          else
            indices.at(D) = 0;
        }
      }
    }

    perversities.emplace_back( aleph::Perversity( values.begin(), values.end() ) );
  }

  return perversities;
}

int main(int argc, char* argv[])
{
  if( argc <= 1 )
    return -1;

  std::string filename = argv[1];
  std::vector<SimplicialComplex> simplicialComplexes;

  aleph::topology::io::LexicographicTriangulationReader reader;
  reader( filename, simplicialComplexes );

  // Create missing faces ----------------------------------------------
  //
  // The triangulations are only specified by their top-level simplices,
  // so they need to be converted before being valid inputs for homology
  // calculations.

  for( auto&& K : simplicialComplexes )
  {
    K.createMissingFaces();
    K.sort();
  }

  // Calculate homology ------------------------------------------------
  //
  // We are only interested in the Betti numbers of the diagrams here as
  // the triangulations are not endowed with any weights or values.

  for( auto&& K : simplicialComplexes )
  {
    bool dualize                    = true;
    bool includeAllUnpairedCreators = true;

    auto diagrams
      = aleph::calculatePersistenceDiagrams( K,
                                            dualize,
                                            includeAllUnpairedCreators );

    for( auto&& D : diagrams )
      std::cout << D.betti() << " ";

    std::cout << "\n";
  }

  // Calculate intersection homology -----------------------------------
  //
  // The basic idea is to first decompose the given simplicial complex
  // into its skeletons. These skeletons then serve as a filtration of
  // the complex. In addition to this, we also calculate a barycentric
  // subdivision of the simplicial complex. The triangulation is hence
  // always "flag-like" following the paper:
  //
  //   Elementary construction of perverse sheaves
  //   Robert MacPhersonl and Kari Vilonen
  //   Inventiones Mathematicae, Volume 84, pp. 403--435, 1986
  //
  // As a last step, we iterate over all possible perversities for the
  // given triangulation and calculate their intersection homology.

  std::vector< std::vector<Signature> > allIntersectionHomologySignatures;
  allIntersectionHomologySignatures.reserve( simplicialComplexes.size() );

  for( auto&& K : simplicialComplexes )
  {
    std::vector<SimplicialComplex> skeletons;
    skeletons.reserve( K.dimension() + 1 );

    aleph::topology::Skeleton skeleton;
    for( std::size_t d = 0; d <= K.dimension(); d++ )
      skeletons.emplace_back( skeleton( d, K ) );

    aleph::topology::BarycentricSubdivision subdivision;
    auto L = subdivision( K );

    // TODO: this is not optimal because the simplicial complexes may
    // share the same dimensionality.
    auto perversities = getPerversities( static_cast<unsigned>( K.dimension() ) );

    std::vector<Signature> signatures;
    signatures.reserve( perversities.size() );

    for( auto&& perversity : perversities )
    {
      auto diagrams  = aleph::calculateIntersectionHomology( L, skeletons, perversity );
      auto signature = makeSignature( diagrams, K.dimension() );

      std::cout << signature << " " << std::flush;

      signatures.push_back( signature );
    }

    std::sort( signatures.begin(), signatures.end() );
    allIntersectionHomologySignatures.emplace_back( signatures );

    std::cout << "\n";
  }
}
