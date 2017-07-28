#include <aleph/topology/io/LexicographicTriangulation.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/PhiPersistence.hh>

#include <aleph/topology/BarycentricSubdivision.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>
#include <aleph/topology/Skeleton.hh>

#include <iostream>
#include <string>
#include <vector>

using DataType          = bool;
using VertexType        = unsigned short;

using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

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

  for( auto&& K : simplicialComplexes )
  {
    std::vector<SimplicialComplex> skeletons;
    skeletons.reserve( K.dimension() + 1 );

    aleph::topology::Skeleton skeleton;
    for( std::size_t d = 0; d <= K.dimension(); d++ )
      skeletons.emplace_back( skeleton( d, K ) );
  }
}
