#include "MainWindow.hh"

#include "persistenceDiagrams/IO.hh"

#include <QAction>
#include <QChart>
#include <QChartView>
#include <QFileDialog>
#include <QLineSeries>
#include <QMdiSubWindow>
#include <QMenuBar>
#include <QScatterSeries>
#include <QStatusBar>

#include <limits>

using namespace QtCharts;

namespace
{

template <class T> QChartView* makeChartView( const aleph::PersistenceDiagram<T>& persistenceDiagram )
{
  QChart* chart                 = new QChart;
  QChartView* chartView         = new QChartView( chart );
  QScatterSeries* scatterSeries = new QScatterSeries;
  QLineSeries* lineSeries       = new QLineSeries;

  T min = std::numeric_limits<T>::max();
  T max = std::numeric_limits<T>::lowest();

  for( auto&& point : persistenceDiagram )
  {
    min = std::min( min, std::min( point.x(), point.y() ) );
    max = std::max( max, std::max( point.x(), point.y() ) );

    scatterSeries->append( point.x(), point.y() );
  }

  // TODO: Make configurable/provide elsewhere
  scatterSeries->setPen( QPen( QColor(196, 30, 58) ) );
  scatterSeries->setBrush( QBrush( QColor( 196, 30, 58 ) ) );
  scatterSeries->setMarkerSize( 1.0 );

  lineSeries->setPen( QPen( Qt::black ) );
  lineSeries->append( min, min );
  lineSeries->append( max, max );

  chart->addSeries( lineSeries );
  chart->addSeries( scatterSeries );
  chart->createDefaultAxes();

  auto legend = chart->legend();
  legend->setVisible( false );

  return chartView;
}

} // anonymous namespace

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
      makeChartView( _persistenceDiagram )
    );

    subWindow->show();
  }
}

} // namespace gui

} // namespace aleph
