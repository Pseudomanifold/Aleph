#include "DataSetItem.hh"
#include "MetaTypes.hh"

#include <QDebug>

namespace aleph
{

namespace gui
{

DataSetItem::DataSetItem( const QString& title,
                          Type type,
                          DataSetItem* parent )
  : _title( title )
  , _parent( parent )
  , _type( type )
{
}

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

  _type = Type::Unspecified;

  if( userType == id_PersistenceDiagram )
  {
    qDebug() << "Identified persistence diagram";
    _type = Type::PersistenceDiagram;
  }
  else if( userType == id_SimplicialComplex )
  {
    qDebug() << "Identified simplicial complex";
    _type = Type::SimplicialComplex;
  }
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
  // TODO: make configurable?
  return 2;
}

QVariant DataSetItem::data( int column ) const
{
  if( column == 0 )
    return _title;

  return QVariant();
}

} // namespace gui

} // namespace aleph
