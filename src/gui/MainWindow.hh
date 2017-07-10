#ifndef ALEPH_GUI_MAIN_WINDOW_HH__
#define ALEPH_GUI_MAIN_WINDOW_HH__

#include <QListView>
#include <QMainWindow>
#include <QMdiArea>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

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
  void createDockWidgets();
  void createMenus();
  void createStatusBar();
  void createToolBars();

private slots:
  void loadPersistenceDiagram();

  void handlePersistenceDiagramClick( const QPointF& point );

private:

  // Widgets -----------------------------------------------------------

  QMdiArea* _mdiArea;
  QListView* _dataSetList;

  // Menus -------------------------------------------------------------

  // Needs to be stored here because different widgets may choose to add
  // themselves to the menu in order to toggle their visibility.
  QMenu* _showMenu;

  // Data --------------------------------------------------------------

  using DataType = double;

  aleph::PersistenceDiagram<DataType> _persistenceDiagram;
};

} // namespace gui

} // namespace aleph

#endif
