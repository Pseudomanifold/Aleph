#include "PersistenceDiagram.hh"

#include <QBrush>
#include <QColor>
#include <QPen>

namespace aleph
{

namespace gui
{

void PersistenceDiagram::setup()
{
  // TODO: Make configurable/provide elsewhere
  _scatterSeries->setPen( QPen( QColor(196, 30, 58) ) );
  _scatterSeries->setBrush( QBrush( QColor( 196, 30, 58 ) ) );
  _scatterSeries->setMarkerSize( 5.0 );

  _lineSeries->setPen( QPen( Qt::black ) );

  auto chart = this->chart();

  chart->addSeries( _lineSeries );
  chart->addSeries( _scatterSeries );
  chart->createDefaultAxes();

  auto legend = chart->legend();
  legend->setVisible( false );

  this->connect( _scatterSeries, SIGNAL( clicked( QPointF ) ),
                 this,           SIGNAL( clicked( QPointF ) ) );
}

} // namespace gui

} // namespace aleph
