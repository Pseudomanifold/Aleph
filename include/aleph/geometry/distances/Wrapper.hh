#ifndef ALEPH_GEOMETRY_DISTANCES_WRAPPER_HH__
#define ALEPH_GEOMETRY_DISTANCES_WRAPPER_HH__

namespace aleph
{

namespace geometry
{

namespace distances
{

/**
  @class Wrapper
  @brief Provides a way to wrap metric calculations

  The purpose of this class is to make it possible to easily wrap metric
  calculations to container classes that have a concept of size. This is
  an easier way of using metrics in different application scenarios. For
  example, this wrapper automatically works with `std::vector`.
*/

template <class Metric, class Container> class Wrapper
{
public:

  Wrapper()
    : _traits( Traits<Metric>() )
  {
  }

  /**
    Main function for calculating the distance between two containers,
    using the specified metric. Note that all trait-based conversions,
    if any, are performed automatically.
  */

  typename Container::value_type operator()( const Container& a, const Container& b )
  {
    return _traits.from( _metric( a.begin(), b.begin(), a.size() ) );
  }

private:
  Metric         _metric; // Metric
  Traits<Metric> _traits; // Traits
};

} // namespace distances

} // namespace geometry

} // namespace aleph

#endif
