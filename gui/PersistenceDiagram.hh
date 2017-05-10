#ifndef ALEPH_GUI_PERSISTENCE_DIAGRAM_HH__
#define ALEPH_GUI_PERSISTENCE_DIAGRAM_HH__

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <QChart>
#include <QChartView>
#include <QLineSeries>
#include <QScatterSeries>

#include <limits>

QT_CHARTS_USE_NAMESPACE

namespace aleph
{

namespace gui
{

class PersistenceDiagram : public QChartView
{

  Q_OBJECT

public:

  template <class T> PersistenceDiagram( const aleph::PersistenceDiagram<T>& persistenceDiagram )
    : QtCharts::QChartView( new QChart )
    , _scatterSeries( new QScatterSeries )
    , _lineSeries( new QLineSeries )
  {
    T min = std::numeric_limits<T>::max();
    T max = std::numeric_limits<T>::lowest();

    for( auto&& point : persistenceDiagram )
    {
      min = std::min( min, std::min( point.x(), point.y() ) );
      max = std::max( max, std::max( point.x(), point.y() ) );

      _scatterSeries->append( point.x(), point.y() );
    }

    _lineSeries->append( min, min );
    _lineSeries->append( max, max );

    this->setup();
  }

signals:
  void clicked( const QPointF );

private:

  void setup();

  QScatterSeries* _scatterSeries; // Points of the persistence diagram
  QLineSeries* _lineSeries;       // Diagonal of the persistence diagram
};

} // namespace gui

} // namespace aleph

#endif
