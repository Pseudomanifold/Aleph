#include <QApplication>

#include "MainWindow.hh"

int main( int argc, char** argv )
{
  QApplication application( argc, argv );

  aleph::gui::MainWindow window;
  window.show();

  return application.exec();
}
