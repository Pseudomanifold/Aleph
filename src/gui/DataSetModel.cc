#include "DataSetModel.hh"
#include "DataSetItem.hh"

#include <QDebug>
#include <QIcon>

namespace aleph
{

namespace gui
{

DataSetModel::DataSetModel( QObject* parent )
  : QAbstractItemModel( parent )
  , _root( new DataSetItem( QString(), QVariant() ) )
{
  QList<DataSetItem*> topLevelItems = {
      new DataSetItem( tr("Persistence diagrams"), QVariant(), _root ),
      new DataSetItem( tr("Point clouds")        , QVariant(), _root ),
      new DataSetItem( tr("Simplicial complexes"), QVariant(), _root )
  };

  foreach( DataSetItem* item, topLevelItems )
  _root->append( item );
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

QVariant DataSetModel::data( const QModelIndex& index, int role ) const
{
  if( !index.isValid() )
    return QVariant();

  if( role == Qt::DisplayRole )
  {
    DataSetItem* item = static_cast<DataSetItem*>( index.internalPointer() );
    return item->data( index.column() );
  }
  else if( role == Qt::DecorationRole && index.column() == 0 )
    return QIcon::fromTheme( "folder" );

  return QVariant();
}

} // namespace gui

} // namespace aleph
