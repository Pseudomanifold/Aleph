#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <aleph/containers/PointCloud.hh>

#include <aleph/config/FLANN.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/RipsExpander.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/math/StepFunction.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/io/SimplicialComplexReader.hh>

#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/PersistenceIndicatorFunction.hh>

#include <aleph/persistenceDiagrams/distances/Bottleneck.hh>
#include <aleph/persistenceDiagrams/distances/Hausdorff.hh>
#include <aleph/persistenceDiagrams/distances/Wasserstein.hh>

#include <aleph/persistenceDiagrams/kernels/MultiScaleKernel.hh>

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/PersistencePairing.hh>

#include <stdexcept>

namespace py = pybind11;

using DataType   = double;
using VertexType = unsigned;

using PointCloud         = aleph::containers::PointCloud<DataType>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using PersistencePairing = aleph::PersistencePairing<VertexType>;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using RipsExpander       = aleph::geometry::RipsExpander<SimplicialComplex>;
using StepFunction       = aleph::math::StepFunction<DataType>;

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
    .def( "__iter__",
      [] ( const PersistenceDiagram& D )
      {
        return py::make_iterator( D.begin(), D.end() );
      }, py::keep_alive<0,1>()
    )
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
    .def_property_readonly( "betti", &PersistenceDiagram::betti )
    .def( "__array__",
      [] (PersistenceDiagram &D)
      {
        auto n_points    = D.size();
        DataType* buffer = new DataType[ 2*n_points ];

        for (auto it = D.begin(); it != D.end(); it++)
        {
          auto index        = std::distance(D.begin(), it);
          buffer[2*index  ] = it->x();
          buffer[2*index+1] = it->y();
        }

        // Based on https://stackoverflow.com/questions/44659924/returning-numpy-arrays-via-pybind11
        // Callback that allows buffer to be freed if not used in python anymore
        py::capsule free_when_done( buffer,
          [] (void* f)
          {
            DataType* buf = reinterpret_cast<DataType*>( f );
            delete[] buf;
          }
        );

        return py::array_t<DataType>(
          {static_cast<unsigned long>(n_points), static_cast<unsigned long>(2)},                          // shape
          {static_cast<unsigned long>(2*sizeof(DataType)), static_cast<unsigned long>(sizeof(DataType))}, // stride
          buffer,                                                                                         // the data pointer
          free_when_done);                                                                                // numpy array references this parent
      }
    );

  using Point = typename PersistenceDiagram::Point;

  py::class_<Point>(m, "PersistenceDiagram.Point" )
    .def( py::init<DataType>() )
    .def( py::init<DataType, DataType>() )
    .def( "__eq__", &Point::operator== )
    .def( "__ne__", &Point::operator!= )
    .def( "__repr__",
      [] ( const Point& p )
      {
        return "<" + std::to_string( p.x() ) + "," + std::to_string( p.y() ) + ">";
      }
    )
    .def_property_readonly( "x", &Point::x )
    .def_property_readonly( "y", &Point::y )
    .def_property_readonly( "persistence", &Point::persistence )
    .def_property_readonly( "unpaired"   , &Point::isUnpaired );
}

void wrapPersistencePairing( py::module& m )
{
  py::class_<PersistencePairing>(m, "PersistencePairing")
    .def( py::init<>() )
    .def( "__bool__",                                   // implicit conversion into booleans
      [] ( const PersistencePairing& pairing )
      {
        return !pairing.empty();
      }
    )
    .def( "__eq__", &PersistencePairing::operator== )   // equality check
    .def( "__ne__", &PersistencePairing::operator!= )   // inequality check
    .def( "__len__", &PersistencePairing::size )        // length
    .def( "__iter__",                                   // iteration
      [] ( const PersistencePairing& pairing )
      {
        return py::make_iterator( pairing.begin(), pairing.end() );
      }, py::keep_alive<0,1>()
    )
    .def( "__repr__",                                   // string-based representation
      [] ( const PersistencePairing& pairing )
      {
        std::ostringstream stream;
        stream << pairing;

        return stream.str();
      }
    );

    //.def( "__array__",
    //  [] (PersistenceDiagram &D) {
    //  auto n_points = D.size();
    //  DataType *buffer = new DataType[2*n_points];

    //  for (auto it = D.begin(); it != D.end(); it++) {
    //    auto index = std::distance(D.begin(), it);
    //    buffer[2*index] = it->x();
    //    buffer[2*index+1] = it->y();
    //  }
    //  // Based on https://stackoverflow.com/questions/44659924/returning-numpy-arrays-via-pybind11
    //  // Callback that allows buffer to be freed if not used in python anymore
    //  py::capsule free_when_done(buffer, [](void *f) {
    //        DataType *buf = reinterpret_cast<DataType *>(f);
    //        delete[] buf;
    //    });

    //  return py::array_t<DataType>(
    //    {static_cast<unsigned long>(n_points), static_cast<unsigned long>(2)}, // shape
    //    {static_cast<unsigned long>(2*sizeof(DataType)), static_cast<unsigned long>(sizeof(DataType))}, // Stride
    //    buffer, // the data pointer
    //    free_when_done); // numpy array references this parent
    //});

  using Point = typename PersistenceDiagram::Point;

  py::class_<Point>(m, "PersistenceDiagram.Point" )
    .def( py::init<DataType>() )
    .def( py::init<DataType, DataType>() )
    .def( "__eq__", &Point::operator== )
    .def( "__ne__", &Point::operator!= )
    .def( "__repr__",
      [] ( const Point& p )
      {
        return "<" + std::to_string( p.x() ) + "," + std::to_string( p.y() ) + ">";
      }
    )
    .def_property_readonly( "x", &Point::x )
    .def_property_readonly( "y", &Point::y )
    .def_property_readonly( "persistence", &Point::persistence )
    .def_property_readonly( "unpaired"   , &Point::isUnpaired );
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

      PointCloud pointCloud(
        static_cast<std::size_t>(n),
        static_cast<std::size_t>(d)
      );

      DataType* target = pointCloud.data();
      DataType* source = reinterpret_cast<DataType*>( bufferInfo.ptr );

      std::copy( source, source + n*d, target );

      using Distance = aleph::geometry::distances::Euclidean<DataType>;
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

void wrapRipsExpander( py::module& m )
{
  py::class_<RipsExpander>(m, "RipsExpander")
    .def( py::init<>() )
    .def( "__call__",
      [] ( RipsExpander& ripsExpander, const SimplicialComplex& K, unsigned dimension )
      {
        return ripsExpander(K, dimension);
      }
    )
    .def( "assignMaximumWeight", &RipsExpander::assignMaximumWeight );
}


void wrapStepFunction( py::module& m )
{
  py::class_<StepFunction>(m, "StepFunction")
    .def( py::init<>() )
    .def( "__init__",
      [] ( StepFunction& instance, const PersistenceDiagram& D )
      {
        new (&instance) StepFunction( aleph::persistenceIndicatorFunction( D ) );
      }
    )
    .def( "__abs__" , &StepFunction::abs )
    .def( "__add__" ,
      [] ( const StepFunction& f, const StepFunction& g )
      {
        return f+g;
      }
    )
    .def( "__sub__" ,
      [] ( const StepFunction& f, const StepFunction& g )
      {
        return f-g;
      }
    )
    .def( "__iadd__", &StepFunction::operator+= )
    .def( "__isub__", &StepFunction::operator-= )
    .def( "__neg__" ,
      [] ( const StepFunction& f )
      {
        return -f;
      }
    )
    .def( "pow",
      [] ( StepFunction& f, double p )
      {
        return f.pow(p);
      }
    )
    .def_property_readonly( "max", &StepFunction::max )
    .def_property_readonly( "sup", &StepFunction::sup )
    .def_property_readonly( "integral", &StepFunction::integral )
    .def( "__call__",
      [] ( const StepFunction& f, DataType x )
      {
        return f(x);
      }
    );

  m.def( "make_persistence_indicator_function",
    [] ( const PersistenceDiagram& diagram )
    {
      return aleph::persistenceIndicatorFunction( diagram );
    }
  );
}

void wrapInputFunctions( py::module& m )
{
  m.def( "load",
    [] ( py::object object )
    {
      std::string filename = py::cast<std::string>( object );

      SimplicialComplex K;

      aleph::topology::io::SimplicialComplexReader reader;
      reader( filename, K);

      return K;
    }
  );

  m.def( "load",
    [] ( py::object object, py::function functor )
    {
      std::string filename = py::cast<std::string>( object );

      SimplicialComplex K;

      aleph::topology::io::SimplicialComplexReader reader;
      reader( filename,
              K,
              [&functor] ( DataType a, DataType b )
              {
                return py::cast<bool>( functor(a,b) );
              }
      );
    }
  );

  m.def( "load_persistence_diagram",
    [] ( py::object object )
    {
      std::string filename = py::cast<std::string>( object );

      return aleph::io::load<DataType>( filename );
    }
  );
}

void wrapNorms( py::module& m )
{
    using namespace pybind11::literals;

    py::module mNorm = m.def_submodule("norms", "Norms on persistence diagrams");
    mNorm.def("totalPersistence",
      [] (const PersistenceDiagram& D, double k, bool weighted)
        {
          return aleph::totalPersistence(D, k, weighted);
        },
      "D"_a,
      "k"_a = 2.0,
      "weighted"_a = false
    );
    mNorm.def("pNorm",
      [] (const PersistenceDiagram& D, double k, bool weighted)
        {
          return aleph::pNorm(D, k, weighted);
        },
      "D"_a,
      "k"_a = 2.0,
      "weighted"_a = false
    );
    mNorm.def("infinityNorm",
      [] (const PersistenceDiagram& D)
        {
          return aleph::infinityNorm(D);
        },
      "D"_a
    );
}


PYBIND11_MODULE(aleph, m)
{
  m.doc() = "Python bindings for Aleph, a library for exploring persistent homology";

  wrapSimplex(m);
  wrapSimplicialComplex(m);
  wrapNorms(m);
  wrapPersistenceDiagram(m);
  wrapPersistencePairing(m);
  wrapPersistentHomologyCalculation(m);
  wrapRipsExpander(m);
  wrapStepFunction(m);
  wrapInputFunctions(m);
}
