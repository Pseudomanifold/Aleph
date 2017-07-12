#include "LoadDataSetHelpers.hh"
#include "MetaTypes.hh"

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <QFileInfo>
#include <QList>

namespace aleph
{

namespace gui
{

QVariant loadData( const QString& file )
{
  QVariant data;

  QList<QString> suffixesPersistenceDiagram = {
    QObject::tr("txt")
  };

  QList<QString> suffixesSimplicialComplex = {
    QObject::tr("gml"),
    QObject::tr("net"),
    QObject::tr("ply"),
    QObject::tr("vtk")
  };

  QFileInfo fileInfo = QFileInfo( file );
  auto suffix        = fileInfo.suffix();

  if( suffixesPersistenceDiagram.contains( suffix ) )
  {
    auto persistenceDiagram = aleph::io::load<DataType>( file.toStdString() );
    data                    = QVariant::fromValue( persistenceDiagram );
  }

  return data;
}

} // namespace gui

} // namespace aleph
