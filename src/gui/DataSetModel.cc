#include "DataSetModel.hh"
#include "DataSetItem.hh"

namespace aleph
{

namespace gui
{

DataSetModel::DataSetModel( QObject* parent )
  : QAbstractItemModel( parent )
  , _root( new DataSetItem( QString(), QVariant() ) )
{
}

DataSetModel::~DataSetModel()
{
  delete _root;
}

QModelIndex DataSetModel::parent( const QModelIndex& index ) const
{
  if( !index.isValid() )
    return QModelIndex();

  auto child  = static_cast<DataSetItem*>( index.internalPointer() );
  auto parent = child->parent();

  if( parent == _root )
    return QModelIndex();

  return this->createIndex( parent->row(), 0, parent );
}

QModelIndex DataSetModel::index( int row, int column,
                                 const QModelIndex& parent ) const
{
  if( !this->hasIndex( row, column, parent ) )
    return QModelIndex();

  DataSetItem* parentItem = nullptr;
  if( !parent.isValid() )
    parentItem = _root;
  else
    parentItem = static_cast<DataSetItem*>( parent.internalPointer() );

  auto child = parentItem->child( row );
  if( child )
    return this->createIndex( row, column, child );
  else
    return QModelIndex();
}

int DataSetModel::rowCount( const QModelIndex& parent ) const
{
  DataSetItem* parentItem = nullptr;
  if( parent.column() > 0 )
    return 0;

  if( !parent.isValid() )
    parentItem = _root;
  else
    parentItem = static_cast<DataSetItem*>( parent.internalPointer() );

  return parentItem->childCount();
}

int DataSetModel::columnCount( const QModelIndex& parent ) const
{
  if( parent.isValid() )
    return static_cast<DataSetItem*>( parent.internalPointer() )->columnCount();
  else
    return _root->columnCount();
}

} // namespace gui

} // namespace aleph
