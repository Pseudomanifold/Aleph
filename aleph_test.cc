#include "BoundaryMatrix.hh"
#include "Vector.hh"
#include "IO.hh"
#include "Standard.hh"
#include "Twist.hh"
#include "PersistencePairs.hh"

#include <iostream>

using namespace aleph;
using namespace representations;

using I  = unsigned;
using V  = Vector<I>;
using BM = BoundaryMatrix<V>;
using SR = StandardReduction;
using TR = TwistReduction;

int main()
{
  auto M = load<BM>( "Triangle.txt" );

  std::cout << "* Boundary matrix\n" << M << "\n"
            << "* Maximum dimension: " << M.getDimension() << "\n";

  computePersistencePairs<SR>( M );
  computePersistencePairs<TR>( M );
}
