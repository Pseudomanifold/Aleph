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

  enum class Type
  {
    PersistenceDiagram,
    SimplicialComplex,
    PointCloud,
    Unspecified
  };

  explicit DataSetItem( const QString& title,
                        Type type,
                        DataSetItem* parent = nullptr );

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
    Returns data stored in the item. The result of the query depends on
    the specified column.
  */

  QVariant data( int column ) const;

private:
  QString  _title;                //< Title of data set (or group)
  QVariant _data;                 //< Data (e.g. persistence diagram)
  QList<DataSetItem*> _children;  //< Children (optional)
  DataSetItem* _parent = nullptr; //< Parent item (optional)

  // Specifies the data type used for the data set item. This is
  // required in order to find the proper item for appending new
  // items into the model.
  Type _type = Type::Unspecified;
};

} // namespace gui

} // namespace aleph

#endif
