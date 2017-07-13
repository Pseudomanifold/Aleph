#ifndef ALEPH_GUI_UNSIGNED_VALIDATOR_HH__
#define ALEPH_GUI_UNSIGNED_VALIDATOR_HH__

#include <QValidator>

namespace aleph
{

namespace gui
{

class UnsignedValidator : public QValidator
{
  Q_OBJECT

public:

  UnsignedValidator( QObject* parent = 0 );

  virtual ~UnsignedValidator()
  {
  }

  /**
    Checks whether the input string is valid. If this is not the case, the
    string will be replaced by an empty string.

    @param string Input string (will be changed)
  */

  virtual void fixup( QString& string ) const;

  /**
    Validates user input. To this end, the function checks whether the
    input string is a valid unsigned value.

    @param input Input string (will not be changed)
    @param pos   Position (unused)

    @returns QValidator::Acceptable or QValidator::Invalid, depending on
    the type check.
  */

  virtual State validate( QString& input , int& pos ) const;
};

} // namespace gui

} // namespace aleph

#endif
