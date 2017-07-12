#include "PersistenceDiagramNormDialog.hh"

#include <QDialogButtonBox>
#include <QDoubleValidator>
#include <QFormLayout>
#include <QRadioButton>
#include <QVBoxLayout>

#include <limits>

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

  int id = 0;
  for( QString s : { tr("Infinity norm"), tr("p-norm"), tr("Total persistence") } )
  {
    QRadioButton* button = new QRadioButton( s );
    button->setChecked( true );

    radioButtonLayout->addWidget( button );

    _normButtonGroup->addButton( button, id++ );
  }

  {
    QDoubleValidator* validator = new QDoubleValidator( this );
    validator->setBottom( 0.0 );

    _powerEdit->setValidator( validator );
  }

  QFormLayout* layout = new QFormLayout;

  layout->addRow( tr("Norm") , radioButtonLayout );
  layout->addRow( tr("Power"), _powerEdit        );

  QDialogButtonBox* buttonBox
    = new QDialogButtonBox(  QDialogButtonBox::Ok
                           | QDialogButtonBox::Cancel );

  QObject::connect( buttonBox, SIGNAL( accepted() ), this, SLOT( accept() ) );
  QObject::connect( buttonBox, SIGNAL( rejected() ), this, SLOT( reject() ) );

  QVBoxLayout* mainLayout = new QVBoxLayout;
  mainLayout->addLayout( layout );
  mainLayout->addWidget( buttonBox );

  this->setLayout( mainLayout );
}

PersistenceDiagramNormDialog::Norm PersistenceDiagramNormDialog::selectedNorm() const
{
  auto selectedButton = _normButtonGroup->checkedId();

  switch( selectedButton )
  {
  case 0:
    return Norm::InfinityNorm;
  case 1:
    return Norm::pNorm;
  case 2:
    return Norm::TotalPersistence;
  }

  return Norm::Undefined;
}

double PersistenceDiagramNormDialog::selectedPower() const
{
  double value = std::numeric_limits<double>::quiet_NaN();
  if( !_powerEdit->text().isEmpty() )
    value = _powerEdit->text().toDouble();

  return value;
}

} // namespace gui

} // namespace aleph
