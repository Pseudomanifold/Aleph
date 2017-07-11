#ifndef ALEPH_GUI_META_TYPES_HH__
#define ALEPH_GUI_META_TYPES_HH__

#include <QMetaType>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

namespace aleph
{

namespace gui
{

using VertexType = unsigned;
using DataType   = double;

using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;

} // namespace gui

} // namespace aleph

Q_DECLARE_METATYPE(aleph::gui::PersistenceDiagram)
Q_DECLARE_METATYPE(aleph::gui::SimplicialComplex)

#endif
