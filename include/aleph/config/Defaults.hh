#ifndef ALEPH_DEFAULTS_HH__
#define ALEPH_DEFAULTS_HH__

#include "persistentHomology/algorithms/Twist.hh"

#include "topology/representations/Vector.hh"

namespace aleph
{

namespace defaults
{

using Index              = unsigned;
using Representation     = topology::representations::Vector<Index>;
using ReductionAlgorithm = persistentHomology::algorithms::Twist;

} // namespace defaults

} // namespace aleph

#endif
