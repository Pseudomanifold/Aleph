#include "SimplicialComplexHelpers.hh"
#include "MetaTypes.hh"

#include <aleph/geometry/RipsExpander.hh>
#include <aleph/geometry/RipsExpanderTopDown.hh>

namespace aleph
{

namespace gui
{

QVariant expandSimplicialComplex( const QVariant& data, unsigned dimension, bool topDown )
{
  QVariant result;

  if( topDown )
  {
    aleph::geometry::RipsExpanderTopDown<SimplicialComplex> expander;
    auto K = data.value<SimplicialComplex>();
    auto L = expander( K, dimension );
    result = QVariant::fromValue( L );

  }
  else
  {
    aleph::geometry::RipsExpander<SimplicialComplex> expander;
    auto K = data.value<SimplicialComplex>();
    auto L = expander( K, dimension );
    result = QVariant::fromValue( L );
  }

  return result;
}

} // namespace gui

} // namespace aleph
