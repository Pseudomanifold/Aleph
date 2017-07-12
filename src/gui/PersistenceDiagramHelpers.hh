#ifndef ALEPH_GUI_PERSISTENCE_DIAGRAM_HELPERS_HH__
#define ALEPH_GUI_PERSISTENCE_DIAGRAM_HELPERS_HH__

#include "PersistenceDiagramNormDialog.hh"

#include <QVariant>

namespace aleph
{

namespace gui
{

double calculateNorm( const QVariant& data,
                      PersistenceDiagramNormDialog::Norm norm,
                      double power );

} // namespace gui

} // namespace aleph

#endif
