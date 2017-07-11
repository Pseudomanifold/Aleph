#ifndef ALEPH_GUI_MAIN_WINDOW_HH__
#define ALEPH_GUI_MAIN_WINDOW_HH__

#include <QMainWindow>
#include <QMdiArea>
#include <QTreeView>

#include "MetaTypes.hh"

namespace aleph
{

namespace gui
{

class DataSetModel;

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

  // Events ------------------------------------------------------------

  void dragEnterEvent( QDragEnterEvent* event );
  void dropEvent( QDropEvent* event );

  // Widgets -----------------------------------------------------------

  QMdiArea* _mdiArea;
  QTreeView* _dataSetView;

  // Menus -------------------------------------------------------------

  // Needs to be stored here because different widgets may choose to add
  // themselves to the menu in order to toggle their visibility.
  QMenu* _showMenu;

  // Data --------------------------------------------------------------

  PersistenceDiagram _persistenceDiagram;

  // Model for storing all data sets that are supplied by the user or
  // created by interactions.
  DataSetModel* _dataSetModel = nullptr;
};

} // namespace gui

} // namespace aleph

#endif
