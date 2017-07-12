#include "MainWindow.hh"

#include "DataSetItem.hh"
#include "DataSetModel.hh"

#include "LoadDataSetHelpers.hh"

#include "PersistenceDiagramHelpers.hh"
#include "PersistenceDiagramNormDialog.hh"
#include "PersistenceDiagramView.hh"

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <QAction>
#include <QDockWidget>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QFileDialog>
#include <QMdiSubWindow>
#include <QMenu>
#include <QMimeData>
#include <QMenuBar>
#include <QStatusBar>

#include <limits>

#include <QDebug>

namespace aleph
{

namespace gui
{

MainWindow::MainWindow()
  : _mdiArea( new QMdiArea( this ) )
  , _dataSetView( new QTreeView( this ) )
  , _dataSetModel( new DataSetModel( this ) )
{
  _dataSetView->setModel( _dataSetModel );
  _dataSetView->setContextMenuPolicy( Qt::CustomContextMenu );

  this->connect( _dataSetView,
                 SIGNAL( customContextMenuRequested(QPoint) ),
                 SLOT( onDataSetContextMenuRequested(QPoint) ) );

  this->createMenus();
  this->createStatusBar();
  this->createToolBars();

  // Needs to be created later on because they modify the menus of the
  // main window.
  this->createDockWidgets();

  // Permits drag & drop events to be handled. I am using this mainly
  // for the quick loading of data sets.
  this->setAcceptDrops( true );

  // MDI area ----------------------------------------------------------

  _mdiArea->setHorizontalScrollBarPolicy( Qt::ScrollBarAsNeeded );
  _mdiArea->setVerticalScrollBarPolicy( Qt::ScrollBarAsNeeded ),

  this->setCentralWidget( _mdiArea );
}

void MainWindow::createDockWidgets()
{
  {
    QDockWidget* dockWidget = new QDockWidget( tr("Data sets"), this );
    dockWidget->setAllowedAreas( Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea );
    dockWidget->setWidget( _dataSetView );

    this->addDockWidget( Qt::LeftDockWidgetArea, dockWidget );

    _showMenu->addAction( dockWidget->toggleViewAction() );
  }

  //
  // FIXME: remove after debugging
  //

  _dataSetModel->add( tr("Iris_dimension_1.txt"), QVariant::fromValue( aleph::io::load<DataType>( "/home/brieck/Projects/Aleph/tests/persistenceDiagrams/Iris_dimension_1.txt" ) ) );
}

void MainWindow::createMenus()
{
  auto menuBar  = this->menuBar();
  auto fileMenu = menuBar->addMenu( tr("&File") );
  _showMenu     = menuBar->addMenu( tr("&Show") );

  // "Load" menu -------------------------------------------------------

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
    _persistenceDiagram = aleph::io::load<DataType>( fileName.toStdString() );

    this->statusBar()->showMessage(
      QString( "Loaded persistence diagram with %1 entries" ).arg( _persistenceDiagram.size() )
    );

    auto pd        = new PersistenceDiagramView( _persistenceDiagram );
    auto subWindow = _mdiArea->addSubWindow( pd );

    this->connect( pd, SIGNAL( clicked( QPointF ) ),
                   this, SLOT( handlePersistenceDiagramClick( QPointF ) ) );

    subWindow->resize( 300, 300 );
    subWindow->show();
  }
}

void MainWindow::handlePersistenceDiagramClick( const QPointF& point )
{
  this->statusBar()->showMessage(
    QString( "Selected point: (%1,%2)" ).arg( point.x() ).arg( point.y() )
  );
}

void MainWindow::onDataSetContextMenuRequested( const QPoint& position )
{
  auto modelIndex = _dataSetView->indexAt( position );

  if( !modelIndex.isValid() || !modelIndex.internalPointer() )
    return;

  DataSetItem* item = static_cast<DataSetItem*>( modelIndex.internalPointer() );
  auto data         = item->data();

  int id_PersistenceDiagram = qMetaTypeId<PersistenceDiagram>();
  int id_SimplicialComplex  = qMetaTypeId<SimplicialComplex>();
  int userType              = data.userType();

  QMenu* contextMenu = new QMenu( this );

  if( userType == id_PersistenceDiagram )
  {
    // Individual persistence diagram selected
    if( item->childCount() == 0 )
    {
      QAction* calculateNormAction    = new QAction( tr("Calculate norm"), contextMenu );
      QAction* visualizeDataSetAction = new QAction( tr("Visualize persistence diagram"), contextMenu );

      this->connect( calculateNormAction   , SIGNAL( triggered() ), SLOT( onCalculateNorm()    ) );
      this->connect( visualizeDataSetAction, SIGNAL( triggered() ), SLOT( onVisualizeDataSet() ) );

      QList<QAction*> actions = {
        calculateNormAction,
        visualizeDataSetAction
      };

      foreach( QAction* action, actions )
      {
        action->setData( data );
        contextMenu->addAction( action );
      }
    }

    // Group of persistence diagrams
    else
    {
      // TODO:
      //  - Comparative visualization
      //  - Distance calculations
    }
  }
  else if( userType == id_SimplicialComplex )
  {
  }

  contextMenu->popup( QCursor::pos() );
}

void MainWindow::onCalculateNorm()
{
  auto action             = qobject_cast<QAction*>( sender() );
  auto data               = action->data();

  PersistenceDiagramNormDialog* dialog = new PersistenceDiagramNormDialog( this );

  auto result = dialog->exec();
  if( result == QDialog::Accepted )
  {
    auto norm   = dialog->selectedNorm();
    auto power  = dialog->selectedPower();
    auto result = calculateNorm( data, norm, power );

    this->statusBar()->showMessage(
      QString( "Persistence diagram norm: %1" ).arg( result )
    );
  }
}

void MainWindow::onVisualizeDataSet()
{
  auto action             = qobject_cast<QAction*>( sender() );
  auto data               = action->data();
  auto persistenceDiagram = data.value<PersistenceDiagram>();
  auto view               = new PersistenceDiagramView( persistenceDiagram );
  auto subWindow          = _mdiArea->addSubWindow( view );

  this->connect( view, SIGNAL( clicked( QPointF ) ),
                 this, SLOT( handlePersistenceDiagramClick( QPointF ) ) );

  subWindow->resize( 300, 300 );
  subWindow->show();
}

void MainWindow::dragEnterEvent( QDragEnterEvent* event )
{
  if( event->mimeData() && event->mimeData()->hasText() )
  {
    event->setDropAction( Qt::CopyAction );
    event->accept();
  }
  else
    event->ignore();
}

void MainWindow::dropEvent( QDropEvent* event )
{
  auto mimeData = event->mimeData();

  if( mimeData->hasUrls() )
  {
    auto urls = mimeData->urls();

    // Ignore an event if multiple URLs are attached to it. We may only handle
    // a single file.
    if( urls.size() > 1 )
    {
      event->ignore();
      return;
    }

    QString file = urls.first().toLocalFile();
    auto data    = loadData( file );

    _dataSetModel->add( file, data );

    this->statusBar()->showMessage( "File: " + file );

    // TODO:
    //  - Check data format (if possible)
    //  - Load it
    //  - Add it to the model
  }
  else
    event->ignore();
}

} // namespace gui

} // namespace aleph
