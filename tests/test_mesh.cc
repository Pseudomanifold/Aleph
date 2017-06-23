#include <tests/Base.hh>

#include <aleph/topology/Mesh.hh>
#include <aleph/topology/MorseSmaleComplex.hh>

#include <vector>

void test1()
{
  ALEPH_TEST_BEGIN( "Simple mesh");

  aleph::topology::Mesh<double> M;

  M.addVertex( 0.0, 0.0, 0.0 );
  M.addVertex( 0.0, 1.0, 0.0 );
  M.addVertex( 1.0, 0.0, 0.0 );
  M.addVertex( 1.5, 1.0, 0.0 );

  std::vector<unsigned> f1 = { 0, 1, 2 };
  std::vector<unsigned> f2 = { 2, 1, 3 };

  M.addFace( f1.begin(), f1.end() );
  M.addFace( f2.begin(), f2.end() );

  ALEPH_ASSERT_EQUAL( M.numVertices(), 4 );
  ALEPH_ASSERT_EQUAL( M.numFaces()   , 2 );

  ALEPH_ASSERT_THROW( M.hasEdge(0,1) );
  ALEPH_ASSERT_THROW( M.hasEdge(1,0) );
  ALEPH_ASSERT_THROW( M.hasEdge(1,2) );
  ALEPH_ASSERT_THROW( M.hasEdge(2,1) );
  ALEPH_ASSERT_THROW( M.hasEdge(0,2) );
  ALEPH_ASSERT_THROW( M.hasEdge(2,0) );

  ALEPH_ASSERT_EQUAL( M.numConnectedComponents(), 1 );

  {
    auto st    = M.star( 0 );
    auto faces = M.faces();

    ALEPH_ASSERT_EQUAL( st.numVertices(), 3 );
    ALEPH_ASSERT_EQUAL( st.numFaces(),    1 );

    ALEPH_ASSERT_THROW( st.hasEdge(0,1) );
    ALEPH_ASSERT_THROW( st.hasEdge(1,0) );
    ALEPH_ASSERT_THROW( st.hasEdge(1,2) );
    ALEPH_ASSERT_THROW( st.hasEdge(2,1) );
    ALEPH_ASSERT_THROW( st.hasEdge(0,2) );
    ALEPH_ASSERT_THROW( st.hasEdge(2,0) );

    ALEPH_ASSERT_EQUAL( faces.front().size(), 3 );

    auto v1 = faces.front().at(0);
    auto v2 = faces.front().at(1);
    auto v3 = faces.front().at(2);

    ALEPH_ASSERT_EQUAL( v1, 0 );
    ALEPH_ASSERT_EQUAL( v2, 1 );
    ALEPH_ASSERT_EQUAL( v3, 2 );
  }

  {
    auto link = M.link( 0 );
    ALEPH_ASSERT_EQUAL( link.size(), 2 );
  }

  ALEPH_TEST_END();
}

void test2()
{
  ALEPH_TEST_BEGIN( "More complex mesh" );

  aleph::topology::Mesh<double> M;

  M.addVertex( 0.0, 0.0, 0.0, 0.0 );
  M.addVertex( 1.0, 0.0, 0.0, 1.0 );
  M.addVertex( 2.0, 0.0, 0.0, 0.0 );
  M.addVertex( 0.0, 1.0, 0.0, 1.0 );
  M.addVertex( 1.0, 1.0, 0.0, 2.0 );
  M.addVertex( 2.0, 1.0, 0.0, 1.0 );
  M.addVertex( 0.0, 2.0, 0.0, 0.0 );
  M.addVertex( 1.0, 2.0, 0.0, 1.0 );
  M.addVertex( 2.0, 2.0, 0.0, 0.0 );

  std::vector<unsigned> f1 = { 0, 1, 4, 3 };
  std::vector<unsigned> f2 = { 1, 2, 5, 4 };
  std::vector<unsigned> f3 = { 4, 5, 8, 7 };
  std::vector<unsigned> f4 = { 3, 4, 7, 6 };

  M.addFace( f1.begin(), f1.end() );
  M.addFace( f2.begin(), f2.end() );
  M.addFace( f3.begin(), f3.end() );
  M.addFace( f4.begin(), f4.end() );

  ALEPH_TEST_END();
}

void test3()
{
  ALEPH_TEST_BEGIN( "Simplicial mesh" );

  aleph::topology::Mesh<double> M;

  M.addVertex( 0.0, 0.0, 0.0, 0.0 );
  M.addVertex( 1.0, 0.0, 0.0, 1.0 );
  M.addVertex( 2.0, 0.0, 0.0, 0.0 );
  M.addVertex( 0.0, 1.0, 0.0, 1.0 );
  M.addVertex( 1.0, 1.0, 0.0, 2.0 );
  M.addVertex( 2.0, 1.0, 0.0, 1.0 );
  M.addVertex( 0.0, 2.0, 0.0, 0.0 );
  M.addVertex( 1.0, 2.0, 0.0, 1.0 );
  M.addVertex( 2.0, 2.0, 0.0, 0.0 );

  std::vector<unsigned> f1 = { 0, 1, 4 };
  std::vector<unsigned> f2 = { 0, 4, 3 };
  std::vector<unsigned> f3 = { 1, 2, 4 };
  std::vector<unsigned> f4 = { 2, 5, 4 };

  std::vector<unsigned> f5 = { 4, 5, 8 };
  std::vector<unsigned> f6 = { 4, 8, 7 };
  std::vector<unsigned> f7 = { 3, 4, 6 };
  std::vector<unsigned> f8 = { 4, 7, 6 };

  M.addFace( f1.begin(), f1.end() );
  M.addFace( f2.begin(), f2.end() );
  M.addFace( f3.begin(), f3.end() );
  M.addFace( f4.begin(), f4.end() );
  M.addFace( f5.begin(), f5.end() );
  M.addFace( f6.begin(), f6.end() );
  M.addFace( f7.begin(), f7.end() );
  M.addFace( f8.begin(), f8.end() );

  {
    auto l1 = M.link(1);
    auto l3 = M.link(3);
    auto l5 = M.link(5);
    auto l7 = M.link(7);

    ALEPH_ASSERT_EQUAL( l1.size(), 3 );
    ALEPH_ASSERT_EQUAL( l1.size(), l3.size() );
    ALEPH_ASSERT_EQUAL( l1.size(), l5.size() );
    ALEPH_ASSERT_EQUAL( l1.size(), l7.size() );

    auto l4 = M.link(4);

    ALEPH_ASSERT_EQUAL( l4.size(), 8 );
  }

  aleph::topology::MorseSmaleComplex<decltype(M)> msc;
  msc( M );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  test1();
  test2();
  test3();
}
