#include "PersistenceDiagram.hh"

#include "distances/Wasserstein.hh"

#include "tests/Base.hh"

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

  {
    auto d11 = wassersteinDistance( D1, D1 );

    assert( d11 >= DataType() );
    assert( d11 == DataType() );

    ALEPH_ASSERT_THROW( d11 >= DataType() );
    ALEPH_ASSERT_THROW( d11 == DataType() );
  }

  Diagram D2;
  D2.add( 0.9, 1.0 );
  D2.add( 1.9, 2.0 );
  D2.add( 2.9, 3.0 );
  D2.add( 3.9, 9.9 );

  {
    auto d12 = wassersteinDistance( D1, D2 );
    auto d21 = wassersteinDistance( D2, D1 );

    ALEPH_ASSERT_THROW( d12 > DataType() );
    ALEPH_ASSERT_THROW( d21 > DataType() );

    ALEPH_ASSERT_THROW( d12 == d21 );
    ALEPH_ASSERT_THROW( std::abs( d12 -  DataType( 3.05 ) ) < 1e-8 );
  }
}
