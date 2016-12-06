#include "algorithms/Standard.hh"
#include "algorithms/Twist.hh"

#include "filtrations/LowerStar.hh"
#include "filtrations/UpperStar.hh"

#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/Calculation.hh"
#include "persistenceDiagrams/Norms.hh"

#include "representations/Set.hh"
#include "representations/Vector.hh"

#include "persistentHomology/Calculation.hh"

#include "topology/BoundaryMatrix.hh"
#include "topology/Conversions.hh"
#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "utilities/String.hh"

#include <iostream>

using namespace aleph;
using namespace algorithms;
using namespace representations;
using namespace filtrations;
using namespace topology;
using namespace utilities;

using I  = unsigned;
using V  = Set<I>;
using BM = BoundaryMatrix<V>;
using SR = Standard;
using TR = Twist;

using S  = Simplex<float, unsigned>;
using SC = SimplicialComplex<S>;

int main()
{
  auto M = BM::load( "Triangle.txt" );

  std::cout << "* Boundary matrix\n" << M << "\n"
            << "* Maximum dimension: " << M.getDimension() << "\n";

  calculatePersistencePairing<SR>( M );
  calculatePersistencePairing<TR>( M );

  calculatePersistencePairing<SR>( M.dualize() );
  calculatePersistencePairing<TR>( M.dualize() );

  std::cout << "* Boundary matrix [doubly-dualized]\n"
            << M.dualize().dualize() << "\n";

  {
    S simplex( {0,1,2} );
    SC K( { {0}, {1}, {2}, {0,1}, {0,2}, {1,2}, {0,1,2} } );

    std::cout << K;

    {
      auto L1 = K;
      auto L2 = K;

      std::vector<float> functionValues
        = { 0.0, 0.0, 1.0, 1.0, 2.0, 3.0, 3.0 };

      LowerStar<S> ls( functionValues.begin(), functionValues.end() );
      UpperStar<S> us( functionValues.begin(), functionValues.end() );

      L1.sort( ls );
      L2.sort( us );

      std::cout << "Lower-star filtration:\n" << L1 << "\n"
                << "Upper-star filtration:\n" << L2 << "\n";
    }

    auto M = makeBoundaryMatrix<BM>( K );

    auto&& pairing1 = calculatePersistencePairing<SR>( M );
    auto&& pairing2 = calculatePersistencePairing<TR>( M );
    auto&& pairing3 = calculatePersistencePairing<SR>( M.dualize() );
    auto&& pairing4 = calculatePersistencePairing<TR>( M.dualize() );

    auto&& diagrams1 = makePersistenceDiagrams( pairing1, K );
    auto&& diagrams2 = makePersistenceDiagrams( pairing2, K );
    auto&& diagrams3 = makePersistenceDiagrams( pairing3, K );
    auto&& diagrams4 = makePersistenceDiagrams( pairing4, K );

    std::cout<< std::string( 80, '-' ) << "\n";
    for( auto&& D : diagrams1 )
      std::cout << D << "\n";
    std::cout<< std::string( 80, '-' ) << "\n";

    for( auto&& D : diagrams2 )
      std::cout << D << "\n";

    std::cout<< std::string( 80, '-' ) << "\n";

    for( auto&& D : diagrams3 )
      std::cout << D << "\n";

    std::cout<< std::string( 80, '-' ) << "\n";

    for( auto&& D : diagrams4 )
      std::cout << D << "\n";

    std::cout<< std::string( 80, '-' ) << "\n";

    for( auto&& D : diagrams1 )
    {
      std::cout << "1-norm:                     " << pNorm( D, 1 ) << "\n"
                << "2-norm:                     " << pNorm( D )    << "\n"
                << "Total degree-1 persistence: " << totalPersistence( D, 1 ) << "\n"
                << "Total degree-2 persistence: " << totalPersistence( D, 2 ) << "\n"
                << "Infinity norm:              " << infinityNorm( D ) << "\n";
    }
  }

  std::string s = " \r\tTest ";
  std::string t = " foo bar   baz\n ";

  std::cout << "#" << trim(s) << "#" << std::endl;

  for( auto&& s : split( t ) )
    std::cout << "*" << s << "*";

  std::cout << std::endl;
}
