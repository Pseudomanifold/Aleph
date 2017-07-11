#include "PersistenceDiagramNormDialog.hh"

#include <QDoubleValidator>
#include <QFormLayout>
#include <QRadioButton>
#include <QVBoxLayout>

namespace aleph
{

namespace gui
{

PersistenceDiagramNormDialog::PersistenceDiagramNormDialog( QWidget* parent )
  : QDialog( parent )
  , _normButtonGroup( new QButtonGroup )
  , _powerEdit( new QLineEdit )
{
  QVBoxLayout* radioButtonLayout = new QVBoxLayout;

  for( QString s : { tr("Infinity norm"), tr("p-norm"), tr("Total persistence") } )
  {
    QRadioButton* button = new QRadioButton( s );
    button->setChecked( true );

    radioButtonLayout->addWidget( button );

    _normButtonGroup->addButton( button );
  }

  {
    QDoubleValidator* validator = new QDoubleValidator( this );
    validator->setBottom( 0.0 );

    _powerEdit->setValidator( validator );
  }

  QFormLayout* layout = new QFormLayout();

  layout->addRow( tr("Norm") , radioButtonLayout );
  layout->addRow( tr("Power"), _powerEdit        );

  this->setLayout( layout );
}

} // namespace gui

} // namespace aleph
