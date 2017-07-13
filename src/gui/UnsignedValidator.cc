#include "UnsignedValidator.hh"

namespace aleph
{

namespace gui
{

UnsignedValidator::UnsignedValidator( QObject* parent )
  : QValidator( parent )
{
}

void UnsignedValidator::fixup( QString& string ) const
{
  int dummy = 0;
  if( this->validate( string, dummy) != QValidator::Acceptable )
    string = "0";
}

QValidator::State UnsignedValidator::validate( QString& input, int& /* pos */ ) const
{
  State result = QValidator::Invalid;

  if( input.length() == 0 )
    result = QValidator::Intermediate;
  else
  {
    bool ok = false;
    input.toUInt( &ok );

    if( ok )
      result = QValidator::Acceptable;
  }

  return result;
}

} // namespace gui

} // namespace aleph
