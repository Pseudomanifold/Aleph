#include <tests/Base.hh>

#include <aleph/geometry/distances/Euclidean.hh>
#include <aleph/geometry/distances/Hamming.hh>
#include <aleph/geometry/distances/Manhattan.hh>

#include <vector>

using namespace aleph;
using namespace geometry;
using namespace distances;

template <class T> void testEuclideanDistance()
{
  ALEPH_TEST_BEGIN( "Euclidean distance" );

  std::vector<T> x = {1,2,3};
  std::vector<T> y = {4,5,6};

  auto functor  = Euclidean<T>();
  auto distance = functor( x.begin(), y.begin(), x.size() );

  ALEPH_ASSERT_THROW( functor.name() == "Euclidean distance" );
  ALEPH_ASSERT_THROW( distance > T(0) );
  ALEPH_ASSERT_EQUAL( distance, T(27) );
  ALEPH_ASSERT_EQUAL( distance, functor( y.begin(), x.begin(), y.size() ) );

  ALEPH_TEST_END();
}

template <class T> void testHammingDistance()
{
  ALEPH_TEST_BEGIN( "Hamming distance" );

  std::vector<T> x = {1,2,3};
  std::vector<T> y = {4,5,6};

  auto functor  = Hamming<T>();
  auto distance = functor( x.begin(), y.begin(), x.size() );

  ALEPH_ASSERT_THROW( functor.name() == "Hamming distance" );
  ALEPH_ASSERT_THROW( distance > T(0) );
  ALEPH_ASSERT_EQUAL( distance, T(3) );
  ALEPH_ASSERT_EQUAL( distance, functor( y.begin(), x.begin(), y.size() ) );

  ALEPH_TEST_END();
}

template <class T> void testManhattanDistance()
{
  ALEPH_TEST_BEGIN( "Manhattan distance" );

  std::vector<T> x = {1,2,3};
  std::vector<T> y = {4,5,6};

  auto functor  = Manhattan<T>();
  auto distance = functor( x.begin(), y.begin(), x.size() );

  ALEPH_ASSERT_THROW( functor.name() == "Manhattan distance" );
  ALEPH_ASSERT_THROW( distance > T(0) );
  ALEPH_ASSERT_EQUAL( distance, T(9) );
  ALEPH_ASSERT_EQUAL( distance, functor( y.begin(), x.begin(), y.size() ) );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  testEuclideanDistance<float> ();
  testEuclideanDistance<double>();

  testHammingDistance<float> ();
  testHammingDistance<double>();

  testManhattanDistance<float> ();
  testManhattanDistance<double>();
}
