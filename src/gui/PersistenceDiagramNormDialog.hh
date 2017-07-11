#ifndef ALEPH_GUI_PERSISTENCE_DIAGRAM_NORM_DIALOG_HH__
#define ALEPH_GUI_PERSISTENCE_DIAGRAM_NORM_DIALOG_HH__

#include <QButtonGroup>
#include <QDialog>
#include <QLineEdit>

namespace aleph
{

namespace gui
{

class PersistenceDiagramNormDialog : public QDialog
{
  Q_OBJECT

public:
  PersistenceDiagramNormDialog( QWidget* parent = nullptr );

private:
  QButtonGroup* _normButtonGroup;
  QLineEdit*    _powerEdit;
};

} // namespace gui

} // namespace aleph

#endif
