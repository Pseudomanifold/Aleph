#ifndef ALEPH_GUI_EXPAND_SIMPLICIAL_COMPLEX_DIALOG_HH__
#define ALEPH_GUI_EXPAND_SIMPLICIAL_COMPLEX_DIALOG_HH__

#include <QDialog>
#include <QLineEdit>

namespace aleph
{

namespace gui
{

class ExpandSimplicialComplexDialog : public QDialog
{
  Q_OBJECT

public:
  ExpandSimplicialComplexDialog( QWidget* parent = nullptr );

  unsigned selectedDimension() const;

private:
  QLineEdit* _dimensionEdit;
};

} // namespace gui

} // namespace aleph

#endif
