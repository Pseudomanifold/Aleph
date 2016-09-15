#ifndef ALEPH_DEFAULTS_HH__
#define ALEPH_DEFAULTS_HH__

#include "algorithms/Twist.hh"
#include "representations/Vector.hh"

namespace aleph
{

namespace defaults
{

using Index              = unsigned;
using Representation     = representations::Vector<Index>;
using ReductionAlgorithm = TwistReduction;

}

}

#endif
