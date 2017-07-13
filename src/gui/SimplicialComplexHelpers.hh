#ifndef ALEPH_GUI_SIMPLICIAL_COMPLEX_HELPERS_HH__
#define ALEPH_GUI_SIMPLICIAL_COMPLEX_HELPERS_HH__

#include <QVariant>

namespace aleph
{

namespace gui
{

QVariant expandSimplicialComplex( const QVariant& data, unsigned dimension, bool topDown = false );

} // namespace gui

} // namespace aleph

#endif
