#include "DataSetItem.hh"
#include "MetaTypes.hh"

#include <QDebug>
#include <QString>

const QList<QString> aleph::gui::DataSetItem::columnNames = {
  QObject::tr("Filename"),
  QObject::tr("Size")
};

namespace aleph
{

namespace gui
{

DataSetItem::DataSetItem( const QString& title,
                          const QVariant& data,
                          DataSetItem* parent )
  : _title( title )
  , _data( data)
  , _parent( parent )
{
  auto id_PersistenceDiagram = qMetaTypeId<PersistenceDiagram>();
  auto id_SimplicialComplex  = qMetaTypeId<SimplicialComplex>();
  auto userType              = data.userType();

  if( userType == id_PersistenceDiagram )
    qDebug() << "Identified persistence diagram";
  else if( userType == id_SimplicialComplex )
    qDebug() << "Identified simplicial complex";
}

DataSetItem::~DataSetItem()
{
  qDeleteAll( _children );
}

void DataSetItem::append( DataSetItem* child )
{
  _children.append( child );
}

int DataSetItem::row() const
{
  if( _parent )
    return _parent->_children.indexOf( const_cast<DataSetItem*>( this ) );

  return 0;
}

int DataSetItem::childCount() const
{
  return _children.count();
}

int DataSetItem::columnCount() const
{
  return DataSetItem::columnNames.count();
}

QVariant DataSetItem::data( int column ) const
{
  switch( column )
  {
  case 0:
    return _title;
  case 1:
    {
      auto id_PersistenceDiagram = qMetaTypeId<PersistenceDiagram>();
      auto id_SimplicialComplex  = qMetaTypeId<SimplicialComplex>();
      auto userType              = _data.userType();

      qulonglong size = 0;

      if( userType == id_PersistenceDiagram )
        size = _data.value<PersistenceDiagram>().size();
      else if( userType == id_SimplicialComplex )
        size = _data.value<SimplicialComplex>().size();

      // Fall back to counting the children of an item if the stored
      // data set is empty
      if( size == 0)
        size = static_cast<qulonglong>( this->childCount() );

      return size;
    }
  }

  return QVariant();
}

} // namespace gui

} // namespace aleph
