#ifndef ALEPH_GUI_DATA_SET_ITEM_HH__
#define ALEPH_GUI_DATA_SET_ITEM_HH__

#include <QList>
#include <QString>
#include <QVariant>

namespace aleph
{

namespace gui
{

class DataSetItem
{
public:

  // Export column names in order to make them usable in the data set
  // model class. I want to keep them close to this class because the
  // indices are used for returning certain data.
  const static QList<QString> columnNames;

  explicit DataSetItem( const QString& title,
                        const QVariant& data,
                        DataSetItem* parent = nullptr );
  ~DataSetItem();

  void append( DataSetItem* child );

  int row() const;

  // Following the Qt nomenclature here, even though I do not really
  // like it.
  int childCount() const;
  int columnCount() const;

  DataSetItem* parent() const         { return _parent; }
  DataSetItem* child( int row ) const { return _children.value( row ); }

  /**
    Returns data type; this corresponds to one of the registered meta
    type IDs.
  */

  int type() const
  {
    return _data.userType();
  }

  /**
    Returns data stored in the item. The result of the query depends on
    the specified column.
  */

  QVariant data( int column ) const;

private:
  QString  _title;                //< Title of data set (or group)
  QVariant _data;                 //< Data (e.g. persistence diagram)
  QList<DataSetItem*> _children;  //< Children (optional)
  DataSetItem* _parent = nullptr; //< Parent item (optional)
};

} // namespace gui

} // namespace aleph

#endif
