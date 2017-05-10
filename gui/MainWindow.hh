#ifndef ALEPH_GUI_MAIN_WINDOW_HH__
#define ALEPH_GUI_MAIN_WINDOW_HH__

#include <QMainWindow>
#include <QMdiArea>

#include "persistenceDiagrams/PersistenceDiagram.hh"

namespace aleph
{

namespace gui
{

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  MainWindow();

private:
  void createMenus();
  void createStatusBar();
  void createToolBars();

private slots:
  void loadPersistenceDiagram();

private:

  // Widgets -----------------------------------------------------------

  QMdiArea* _mdiArea;

  // Data --------------------------------------------------------------

  using DataType = double;

  aleph::PersistenceDiagram<DataType> _persistenceDiagram;
};

} // namespace gui

} // namespace aleph

#endif
