#ifndef ALEPH_GUI_LOAD_DATA_SET_HELPERS_HH__
#define ALEPH_GUI_LOAD_DATA_SET_HELPERS_HH__

#include <QString>
#include <QVariant>

namespace aleph
{

namespace gui
{

QVariant loadData( const QString& file );
QVariant loadPersistenceDiagram( const QString& file );

} // namespace gui

} // namespace aleph

#endif
