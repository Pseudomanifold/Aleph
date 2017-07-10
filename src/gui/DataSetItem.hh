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
  explicit DataSetItem( const QString& title,
                        const QVariant& data,
                        DataSetItem* parent = nullptr );
  ~DataSetItem();

  void append( DataSetItem* child );

  int children() const;
  int row() const;

  DataSetItem* parent() const         { return _parent; }
  DataSetItem* child( int row ) const { return _children.value( row ); }

  QString title() const               { return _title; }
  QVariant data() const               { return _data;  }

private:
  QString  _title;                //< Title of data set (or group)
  QVariant _data;                 //< Data (e.g. persistence diagram)
  QList<DataSetItem*> _children;  //< Children (optional)

  DataSetItem* _parent = nullptr; //< Parent item (optional)
};

} // namespace gui

} // namespace aleph

#endif
