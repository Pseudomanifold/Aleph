#include "ExpandSimplicialComplexDialog.hh"
#include "UnsignedValidator.hh"

#include <QDialogButtonBox>
#include <QLabel>
#include <QVBoxLayout>

namespace aleph
{

namespace gui
{

ExpandSimplicialComplexDialog::ExpandSimplicialComplexDialog( QWidget* parent )
  : QDialog( parent )
  , _dimensionEdit( new QLineEdit( this ) )
{
  _dimensionEdit->setValidator( new UnsignedValidator );
  _dimensionEdit->setPlaceholderText( tr("Dimension") );

  QLabel* label
    = new QLabel( tr("Please enter the maximum dimension for expanding the simplicial\n")
                 +tr("complex. A value of <em>d</em> ensures that the complex contains\n")
                 +tr("<em>d</em>-simplices if possible.") );

  label->setTextFormat( Qt::RichText );
  label->setWordWrap( true );

  QDialogButtonBox* buttonBox
    = new QDialogButtonBox(  QDialogButtonBox::Ok
                           | QDialogButtonBox::Cancel );

  QObject::connect( buttonBox, SIGNAL( accepted() ), this, SLOT( accept() ) );
  QObject::connect( buttonBox, SIGNAL( rejected() ), this, SLOT( reject() ) );

  QVBoxLayout* mainLayout = new QVBoxLayout( this );
  mainLayout->addWidget( label );
  mainLayout->addWidget( _dimensionEdit );
  mainLayout->addWidget( buttonBox );

  this->setLayout( mainLayout );
}

unsigned ExpandSimplicialComplexDialog::selectedDimension() const
{
  auto text = _dimensionEdit->text();

  if( text.isEmpty() )
    return 0;
  else
    return text.toUInt();
}

} // namespace gui

} // namespace aleph
