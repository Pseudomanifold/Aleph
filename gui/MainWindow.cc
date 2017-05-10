#include "gui/MainWindow.hh"
#include "gui/PersistenceDiagram.hh"

#include "persistenceDiagrams/IO.hh"

#include <QAction>
#include <QFileDialog>
#include <QMdiSubWindow>
#include <QMenuBar>
#include <QStatusBar>

#include <limits>

namespace aleph
{

namespace gui
{

MainWindow::MainWindow()
  : _mdiArea( new QMdiArea( this ) )
{
  this->createMenus();
  this->createStatusBar();
  this->createToolBars();

  // MDI area ----------------------------------------------------------

  _mdiArea->setHorizontalScrollBarPolicy( Qt::ScrollBarAsNeeded );
  _mdiArea->setVerticalScrollBarPolicy( Qt::ScrollBarAsNeeded ),

  this->setCentralWidget( _mdiArea );
}

void MainWindow::createMenus()
{
  auto menuBar  = this->menuBar();
  auto fileMenu = menuBar->addMenu( tr("&File") );
  auto loadMenu = fileMenu->addMenu( tr("Load") );

  QAction* loadPersistenceDiagram
    = new QAction( "Persistence diagram", loadMenu );

  loadMenu->addAction( loadPersistenceDiagram );

  this->connect( loadPersistenceDiagram, &QAction::triggered, this, &MainWindow::loadPersistenceDiagram );
}

void MainWindow::createStatusBar()
{
  this->statusBar()->showMessage( tr("Welcome!"), 2000 );
}

void MainWindow::createToolBars()
{
}

void MainWindow::loadPersistenceDiagram()
{
  auto fileName = QFileDialog::getOpenFileName( this );
  if( !fileName.isEmpty() )
  {
    _persistenceDiagram = aleph::load<DataType>( fileName.toStdString() );

    this->statusBar()->showMessage(
      QString( "Loaded persistence diagram with %1 entries" ).arg( _persistenceDiagram.size() )
    );

    auto subWindow = _mdiArea->addSubWindow(
      new PersistenceDiagram( _persistenceDiagram )
    );

    subWindow->show();
  }
}

} // namespace gui

} // namespace aleph
