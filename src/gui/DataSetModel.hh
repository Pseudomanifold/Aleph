#ifndef ALEPH_GUI_DATA_SET_MODEL_HH__
#define ALEPH_GUI_DATA_SET_MODEL_HH__

#include <QAbstractItemModel>
#include <QObject>

#if 0
class TreeModel : public QAbstractItemModel
{
    Q_OBJECT

public:
    explicit TreeModel(const QString &data, QObject *parent = 0);
    ~TreeModel();

    QVariant data(const QModelIndex &index, int role) const override;
    Qt::ItemFlags flags(const QModelIndex &index) const override;
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const override;
    QModelIndex index(int row, int column,
                      const QModelIndex &parent = QModelIndex()) const override;
    QModelIndex parent(const QModelIndex &index) const override;
    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;

private:
    void setupModelData(const QStringList &lines, TreeItem *parent);

    TreeItem *rootItem;
};
#endif

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

private:

  // This is the root of all data items. It does *not* contain any
  // additional data and is generally not shown directly.
  DataSetItem* _root = nullptr;
};

} // namespace gui

} // namespace aleph

#endif
