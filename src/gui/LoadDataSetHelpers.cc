#include "LoadDataSetHelpers.hh"
#include "MetaTypes.hh"

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <aleph/topology/io/GML.hh>
#include <aleph/topology/io/Pajek.hh>
#include <aleph/topology/io/PLY.hh>
#include <aleph/topology/io/VTK.hh>

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
  else if( suffixesSimplicialComplex.contains( suffix ) )
  {
    SimplicialComplex K;
    if( suffix == QObject::tr("gml") )
    {
      aleph::topology::io::GMLReader reader;
      reader( file.toStdString(), K );

    }
    else if( suffix == QObject::tr("net") )
    {
      aleph::topology::io::PajekReader reader;
      reader( file.toStdString(), K );
    }
    else if( suffix == QObject::tr("ply") )
    {
      aleph::topology::io::PLYReader reader;
      reader( file.toStdString(), K );
    }
    else if( suffix == QObject::tr("vtk") )
    {
      aleph::topology::io::VTKStructuredGridReader reader;
      reader( file.toStdString(), K );
    }

   data = QVariant::fromValue( K );
  }

  return data;
}

} // namespace gui

} // namespace aleph
