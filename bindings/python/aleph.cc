#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <aleph/containers/PointCloud.hh>

#include <aleph/config/FLANN.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/MultiScaleKernel.hh>

#include <aleph/persistenceDiagrams/distances/Bottleneck.hh>
#include <aleph/persistenceDiagrams/distances/Hausdorff.hh>
#include <aleph/persistenceDiagrams/distances/Wasserstein.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <stdexcept>

namespace py = pybind11;

using DataType   = double;
using VertexType = unsigned;

using PointCloud         = aleph::containers::PointCloud<DataType>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;

#ifdef ALEPH_WITH_FLANN
  template <class Distance> using NearestNeighbours = aleph::geometry::FLANN<PointCloud, Distance>;

#else
  template <class Distance> using NearestNeighbours = aleph::geometry::BruteForce<PointCloud, Distance>;
#endif

void wrapSimplex( py::module& m )
{
  py::class_<Simplex>(m, "Simplex")
    .def( py::init<>() )
    .def( py::init<VertexType, DataType>() )
    .def( py::init<Simplex, DataType>() )
    .def( "__init__",
          [] ( Simplex& instance, py::list vertices_ )
          {
            std::vector<VertexType> vertices;
            for( auto handle : vertices_ )
              vertices.push_back( py::cast<VertexType>( handle ) );

            new (&instance) Simplex( vertices.begin(), vertices.end() );
          }
    )
    .def( "__init__",
          [] ( Simplex& instance, py::list vertices_, DataType dataType )
          {
            std::vector<VertexType> vertices;
            for( auto handle : vertices_ )
              vertices.push_back( py::cast<VertexType>( handle ) );

            new (&instance) Simplex( vertices.begin(), vertices.end(), dataType );
          }
    )
    .def( "__bool__",
          [] ( const Simplex& simplex )
          {
            return !simplex.empty();
          }
    )
    .def( "__contains__",
          [] ( const Simplex& simplex, VertexType v )
          {
            return simplex.contains(v);
          }
    )
    .def( "__getitem__",
          [] ( const Simplex& simplex, VertexType i )
          {
            return simplex[i];
          }
    )
    .def( "__iter__",
          [] ( const Simplex& simplex )
          {
            return py::make_iterator( simplex.begin(), simplex.end() );
          }, py::keep_alive<0,1>()
    )
    .def( "__eq__" , &Simplex::operator== )
    .def( "__ne__" , &Simplex::operator!= )
    .def( "__lt__" , &Simplex::operator< )
    .def( "__len__", &Simplex::size )
    .def( "__repr__",
          [] ( const Simplex& simplex )
          {
            std::ostringstream stream;
            stream << simplex;

            return stream.str();
          }
    )
    .def( "__reversed__",
          [] ( const Simplex& simplex )
          {
            return py::make_iterator( simplex.rbegin(), simplex.rend() );
          }, py::keep_alive<0,1>()
    )
    .def_property_readonly( "boundary",
          [] ( const Simplex& simplex )
          {
            return py::make_iterator( simplex.begin_boundary(), simplex.end_boundary() );
          }, py::keep_alive<0,1>()
    )
    .def_property_readonly( "dimension", &Simplex::dimension )
    .def_property_readonly( "empty"    , &Simplex::empty )
    .def_property("data"  , &Simplex::data, &Simplex::setData )
    .def_property("weight", &Simplex::data, &Simplex::setData );
}

void wrapSimplicialComplex( py::module& m )
{
  py::class_<SimplicialComplex>(m, "SimplicialComplex")
    .def( py::init<>() )
    .def( "__init__",
      [] ( SimplicialComplex& instance, py::list simplices_ )
      {
        std::vector<Simplex> simplices;
        for( auto simplexHandle : simplices_ )
        {
          // Let us first try to obtain a simplex from each handle
          // in order to rapidly build a complex.
          try
          {
            auto simplex = py::cast<Simplex>( simplexHandle );
            simplices.push_back( simplex );
          }

          // Assume that the list contains only lists of vertices
          // and convert them directly to simplices.
          catch( py::cast_error& )
          {
            std::vector<VertexType> vertices;
            DataType data = DataType();

            try
            {
              auto&& vertices_ = py::cast<py::list>( simplexHandle );

              for( auto vertexHandle : vertices_ )
                vertices.push_back( py::cast<VertexType>( vertexHandle ) );
            }
            catch( py::cast_error& )
            {
              auto&& tuple_    = py::cast<py::tuple>( simplexHandle );

              if( tuple_.size() != 2 )
                throw std::runtime_error( "Unsupported number of tuple elements" );

              auto&& vertices_ = py::cast<py::list>( tuple_[0] );
              data             = py::cast<DataType>( tuple_[1] );

              for( auto vertexHandle : vertices_ )
                vertices.push_back( py::cast<VertexType>( vertexHandle ) );
            }

            simplices.push_back( Simplex( vertices.begin(), vertices.end(), data ) );
          }

        }

        new (&instance) SimplicialComplex( simplices.begin(), simplices.end() );
      }
    )
    .def( "__bool__",
      [] ( const SimplicialComplex& K )
      {
        return !K.empty();
      }
    )
    .def( "__contains__", &SimplicialComplex::contains )
    .def( "__getitem__",  &SimplicialComplex::operator[] ) // FIXME: might want to handle negative indices?
    .def( "__iter__",
      [] ( const SimplicialComplex& K )
      {
        return py::make_iterator( K.begin(), K.end() );
      }, py::keep_alive<0,1>()
    )
    .def( "__len__", &SimplicialComplex::size )
    .def( "__repr__",
      [] ( const SimplicialComplex& K )
      {
        std::ostringstream stream;
        stream << K;

        return stream.str();
      }
    )
    .def( "append", &SimplicialComplex::push_back )
    .def( "append",
      [] ( SimplicialComplex& K, py::list vertices_ )
      {
        std::vector<VertexType> vertices;
        for( auto vertexHandle : vertices_ )
          vertices.push_back( py::cast<VertexType>( vertexHandle ) );

        K.push_back( Simplex( vertices.begin(), vertices.end() ) );
      }
    )
    .def( "sort",
      [] ( SimplicialComplex& K )
      {
        K.sort();
        return K;
      }
    )
    .def( "sort",
      [] ( SimplicialComplex& K, py::function functor )
      {
        K.sort(
          [&functor] ( const Simplex& s, const Simplex& t )
          {
            return py::cast<bool>( functor(s,t) );
          }
        );

        return K;
      }
    )
    .def_property_readonly( "dimension", &SimplicialComplex::dimension );
}

void wrapPersistenceDiagram( py::module& m )
{
  py::class_<PersistenceDiagram>(m, "PersistenceDiagram")
    .def( py::init<>() )
    .def( "__bool__",
      [] ( const PersistenceDiagram& D )
      {
        return !D.empty();
      }
    )
    .def( "__eq__", &PersistenceDiagram::operator== )
    .def( "__ne__", &PersistenceDiagram::operator!= )
    .def( "__len__", &PersistenceDiagram::size )
    .def( "__repr__",
      [] ( const PersistenceDiagram& D )
      {
        std::ostringstream stream;
        stream << D;

        return stream.str();
      }
    )
    .def( "removeDiagonal", &PersistenceDiagram::removeDiagonal )
    .def( "removeUnpaired", &PersistenceDiagram::removeUnpaired )
    .def_property( "dimension", &PersistenceDiagram::setDimension, &PersistenceDiagram::dimension )
    .def_property_readonly( "betti", &PersistenceDiagram::betti );
}

void wrapPersistentHomologyCalculation( py::module& m )
{
  using namespace pybind11::literals;

  m.def( "calculatePersistenceDiagrams",
    [] ( const SimplicialComplex& K )
    {
      return aleph::calculatePersistenceDiagrams( K );
    }
  );

  m.def( "calculatePersistenceDiagrams",
    [] ( py::buffer buffer, DataType epsilon, unsigned dimension )
    {
      py::buffer_info bufferInfo = buffer.request();

      if( bufferInfo.ndim != 2 || bufferInfo.shape.size() != 2 )
        throw std::runtime_error( "Only two-dimensional buffers are supported" );

      if( bufferInfo.format != py::format_descriptor<DataType>::format() )
        throw std::runtime_error( "Unexpected format" );

      auto n = bufferInfo.shape[0];
      auto d = bufferInfo.shape[1];

      PointCloud pointCloud(n,d);

      DataType* target = pointCloud.data();
      DataType* source = reinterpret_cast<DataType*>( bufferInfo.ptr );

      std::copy( source, source + n*d, target );

      using Distance = aleph::distances::Euclidean<DataType>;
      dimension      = dimension > 0 ? dimension : static_cast<unsigned>( pointCloud.dimension() + 1 );

      auto K         = aleph::geometry::buildVietorisRipsComplex(
        NearestNeighbours<Distance>( pointCloud ),
        epsilon,
        dimension
      );

      return aleph::calculatePersistenceDiagrams( K );
    },
    "buffer"_a,
    "epsilon"_a   = DataType(),
    "dimension"_a = 0
  );
}

void wrapDistanceCalculations( py::module& m )
{
  using namespace pybind11::literals;

  m.def( "bottleneckDistance",
    [] (const PersistenceDiagram& D1, const PersistenceDiagram& D2 )
    {
      return aleph::distances::bottleneckDistance( D1, D2 );
    }
  );

  m.def( "hausdorffDistances",
    [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2 )
    {
      return aleph::distances::hausdorffDistance( D1, D2 );
    }
  );

  m.def( "wassersteinDistance",
    [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, DataType p )
    {
      return aleph::distances::wassersteinDistance( D1, D2, p );
    },
    "D1"_a,
    "D2"_a,
    "p"_a = DataType(1)
  );
}

void wrapKernelCalculations( py::module& m )
{
  using namespace pybind11::literals;

  m.def( "multiScaleKernel",
    [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, double sigma )
    {
      return aleph::multiScaleKernel( D1, D2, sigma );
    }
  );

  m.def( "multiScalePseudoMetric",
    [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, double sigma )
    {
      return aleph::multiScalePseudoMetric( D1, D2, sigma );
    }
  );
}

PYBIND11_PLUGIN(aleph)
{
  py::module m("aleph", "Python bindings for Aleph, a library for exploring persistent homology");

  wrapSimplex(m);
  wrapSimplicialComplex(m);
  wrapPersistenceDiagram(m);
  wrapPersistentHomologyCalculation(m);

  return m.ptr();
}
