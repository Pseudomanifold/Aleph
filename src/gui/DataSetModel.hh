#ifndef ALEPH_GUI_DATA_SET_MODEL_HH__
#define ALEPH_GUI_DATA_SET_MODEL_HH__

#include <QAbstractItemModel>
#include <QObject>

namespace aleph
{

namespace gui
{

class DataSetItem;

class DataSetModel : public QAbstractItemModel
{
  Q_OBJECT

public:
  DataSetModel( QObject* parent = nullptr );
  ~DataSetModel();

  QModelIndex parent( const QModelIndex& index ) const override;
  QModelIndex index( int row, int column, const QModelIndex& parent = QModelIndex() ) const override;

  int rowCount( const QModelIndex& parent = QModelIndex() ) const override;
  int columnCount( const QModelIndex& parent = QModelIndex() ) const override;

  QVariant data( const QModelIndex& index, int role ) const override;

  /**
    Main entry function for adding new data sets to the model. The model
    will automatically check the type of the data set and insert it
    accordingly.

    Unknown data sets will be refused.
  */

  void add( const QString& title, const QVariant& data );

private:

  // This is the root of all data items. It does *not* contain any
  // additional data and is generally not shown directly.
  DataSetItem* _root = nullptr;
};

} // namespace gui

} // namespace aleph

#endif
