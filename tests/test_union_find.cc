#include <tests/Base.hh>

#include <aleph/topology/UnionFind.hh>

#include <iterator>
#include <set>
#include <typeinfo>
#include <vector>

using namespace aleph::topology;

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Union--Find (" + std::string( typeid(T).name() ) + ")" );

  std::vector<T> vertices = {1,2,3,4,5,6,7,8};

  UnionFind<T> uf( vertices.begin(), vertices.end() );

  for( auto&& vertex : vertices )
    ALEPH_ASSERT_EQUAL( uf.find(vertex), vertex );

  uf.merge(1,2);
  uf.merge(5,6);
  uf.merge(5,8);

  ALEPH_ASSERT_EQUAL( uf.find(1), uf.find(2) );

  ALEPH_ASSERT_EQUAL( uf.find(5), uf.find(6) );
  ALEPH_ASSERT_EQUAL( uf.find(6), uf.find(5) );
  ALEPH_ASSERT_EQUAL( uf.find(5), uf.find(8) );
  ALEPH_ASSERT_EQUAL( uf.find(8), uf.find(5) );

  uf.merge(3,4);
  uf.merge(1,5);

  ALEPH_ASSERT_EQUAL( uf.find(3), uf.find(4) );
  ALEPH_ASSERT_EQUAL( uf.find(7), 7          );

  std::set<T> roots;
  uf.roots( std::inserter( roots, roots.begin() ) );

  ALEPH_ASSERT_EQUAL( roots.size(), 3 );
  ALEPH_ASSERT_THROW( roots.find(4) != roots.end() );
  ALEPH_ASSERT_THROW( roots.find(7) != roots.end() );
  ALEPH_ASSERT_THROW( roots.find(8) != roots.end() );

  std::set<T> component1;
  std::set<T> component2;
  std::set<T> component3;

  uf.get( 4, std::inserter( component1, component1.begin() ) );
  uf.get( 7, std::inserter( component2, component2.begin() ) );
  uf.get( 8, std::inserter( component3, component3.begin() ) );

  ALEPH_ASSERT_EQUAL( component1.size(), 2 );
  ALEPH_ASSERT_EQUAL( component2.size(), 1 );
  ALEPH_ASSERT_EQUAL( component3.size(), 5 );

  ALEPH_ASSERT_THROW( component1 == std::set<T>( {3,4}       ) );
  ALEPH_ASSERT_THROW( component2 == std::set<T>( {7}         ) );
  ALEPH_ASSERT_THROW( component3 == std::set<T>( {1,2,5,6,8} ) );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  test<unsigned short>();
  test<short>         ();
  test<int>           ();
  test<unsigned>      ();
  test<long>          ();
  test<unsigned long> ();
}
