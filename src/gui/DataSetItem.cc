#include "DataSetItem.hh"

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

} // namespace gui

} // namespace aleph
