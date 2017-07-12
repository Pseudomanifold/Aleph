#include "PersistenceDiagramHelpers.hh"
#include "MetaTypes.hh"

#include <aleph/persistenceDiagrams/Norms.hh>

#include <limits>

namespace aleph
{

namespace gui
{

double calculateNorm( const QVariant& data,
                      PersistenceDiagramNormDialog::Norm norm,
                      double power )
{
  using Norm              = PersistenceDiagramNormDialog::Norm;
  auto persistenceDiagram = data.value<PersistenceDiagram>();
  double result           = std::numeric_limits<double>::quiet_NaN();

  switch( norm )
  {
    case Norm::InfinityNorm:
      result = static_cast<double>( aleph::infinityNorm( persistenceDiagram ) );
      break;

    case Norm::pNorm:
      result = static_cast<double>( aleph::pNorm( persistenceDiagram, power ) );
      break;

    case Norm::TotalPersistence:
      result = static_cast<double>( aleph::totalPersistence( persistenceDiagram, power ) );
      break;

    default:
      break;
  }

  return result;
}

} // namespace gui

} // namespace aleph
