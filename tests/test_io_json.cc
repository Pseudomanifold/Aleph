#include <tests/Base.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/io/JSON.hh>

#include <fstream>
#include <vector>

template <class T> void testReading()
{
  auto filename            = CMAKE_SOURCE_DIR + std::string("/tests/input/Persistence_diagrams.json");
  auto diagrams            = aleph::io::readJSON<T>( filename );

  ALEPH_ASSERT_EQUAL( diagrams.size()    ,  5 );
  ALEPH_ASSERT_EQUAL( diagrams[0].size() , 38 );
  ALEPH_ASSERT_EQUAL( diagrams[0].betti(), 12 );
  ALEPH_ASSERT_EQUAL( diagrams[3].size(),   0 );
  ALEPH_ASSERT_EQUAL( diagrams[4].size(),   0 );
}

int main(int, char**)
{
  testReading<double>();
  testReading<float >();
}
