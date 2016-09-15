#include "BoundaryMatrix.hh"
#include "IO.hh"
#include "PersistencePairs.hh"

#include "Dualization.hh"
#include "PersistenceDiagram.hh"
#include "PersistenceDiagramCalculation.hh"
#include "PersistenceDiagramNorms.hh"
#include "Simplex.hh"
#include "SimplicialComplex.hh"
#include "SimplicialComplexConversions.hh"

#include "algorithms/Standard.hh"
#include "algorithms/Twist.hh"

#include "filtrations/LowerStar.hh"
#include "filtrations/UpperStar.hh"

#include "utilities/String.hh"

#include "representations/Vector.hh"

#include <iostream>

using namespace aleph;
using namespace representations;
using namespace filtrations;
using namespace utilities;

using I  = unsigned;
using V  = Vector<I>;
using BM = BoundaryMatrix<V>;
using SR = StandardReduction;
using TR = TwistReduction;

using S  = Simplex<float, unsigned>;
using SC = SimplicialComplex<S>;

int main()
{
  auto M = load<BM>( "Triangle.txt" );

  std::cout << "* Boundary matrix\n" << M << "\n"
            << "* Maximum dimension: " << M.getDimension() << "\n";

  computePersistencePairs<SR>( M );
  computePersistencePairs<TR>( M );

  computePersistencePairs<SR>( dualize( M ) );
  computePersistencePairs<TR>( dualize( M ) );

  std::cout << "* Boundary matrix [doubly-dualized]\n"
            << dualize( dualize( M ) ) << "\n";

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

    auto&& pairing1 = computePersistencePairs<SR>( M );
    auto&& pairing2 = computePersistencePairs<TR>( M );
    auto&& pairing3 = computePersistencePairs<SR>( dualize( M ) );
    auto&& pairing4 = computePersistencePairs<TR>( dualize( M ) );

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
