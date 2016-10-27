#include <iostream>

#include "PersistenceDiagram.hh"

#include "distances/Wasserstein.hh"

using namespace aleph;
using namespace aleph::distances;

using DataType = double;
using Diagram  = PersistenceDiagram<DataType>;

int main()
{
  Diagram D1;
  D1.add( 0.9, 1.0 );
  D1.add( 1.9, 2.0 );
  D1.add( 2.9, 3.0 );
  D1.add( 3.9, 4.0 );

  Diagram D2;
  D2.add( 0.9, 1.0 );
  D2.add( 1.9, 2.0 );
  D2.add( 2.9, 3.0 );
  D2.add( 3.9, 4.0 );

  auto d11 = wassersteinDistance( D1, D1 );
  auto d12 = wassersteinDistance( D1, D2 );

  std::cerr << "d11 = " << d11 << std::endl
            << "d12 = " << d12 << std::endl;
}
