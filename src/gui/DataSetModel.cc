#include "DataSetModel.hh"

#include "DataSetItem.hh"
#include "MetaTypes.hh"

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
      new DataSetItem( tr("Persistence diagrams"), QVariant::fromValue( PersistenceDiagram() ), _root ),
      new DataSetItem( tr("Point clouds")        , QVariant::fromValue( SimplicialComplex()  ), _root ),
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

QVariant DataSetModel::headerData( int section,
                                   Qt::Orientation orientation,
                                   int role ) const
{
  if( orientation == Qt::Horizontal && role == Qt::DisplayRole )
    return DataSetItem::columnNames[section];

  return QVariant();
}

QVariant DataSetModel::data( const QModelIndex& index, int role ) const
{
  if( !index.isValid() )
    return QVariant();

  DataSetItem* item = static_cast<DataSetItem*>( index.internalPointer() );

  if( role == Qt::DecorationRole && index.column() == 0 && item->parent() == _root )
    return QIcon::fromTheme( "folder" );
  else if( role == Qt::DisplayRole )
    return item->data( index.column() );

  return QVariant();
}

void DataSetModel::add( const QString& title, const QVariant& data )
{
  int id_PersistenceDiagram = qMetaTypeId<PersistenceDiagram>();
  int id_SimplicialComplex  = qMetaTypeId<SimplicialComplex>();
  int userType              = data.userType();

  if( userType != id_PersistenceDiagram && userType != id_SimplicialComplex )
  {
    qDebug() << "Ignoring unknown user type";
    return;
  }

  for( int i = 0; i < _root->childCount(); i++ )
  {
    auto child = _root->child( i );
    if( child && child->type() == userType )
    {
      qDebug() << "Found proper parent item; adding data";
      child->append( new DataSetItem( title, data, child ) );
    }
  }
}

} // namespace gui

} // namespace aleph
