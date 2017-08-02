#include <tests/Base.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/distances/Euclidean.hh>
#include <aleph/geometry/WitnessComplex.hh>

template <class T> test()
{
  using PointCloud = aleph::PointCloud<T>;
}

int main(int, char**)
{
  test<float> ();
  test<double>();
}
