/**
  @file  aleph.cc
  @brief Main file for Aleph's Python bindings
*/

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Configuration -------------------------------------------------------

#include <aleph/config/Defaults.hh>
#include <aleph/config/FLANN.hh>

// Containers ----------------------------------------------------------

#include <aleph/containers/PointCloud.hh>

// Geometry ------------------------------------------------------------
//
// Various classes for nearest-neighbour calculations and the expansion
// process of a simplicial complex.

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/RipsExpander.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

// Step functions ------------------------------------------------------

#include <aleph/math/StepFunction.hh>

// Simplicial complexes ------------------------------------------------
//
// Main simplex class, simplicial complexes, filtrations, and the
// assorted reader classes.

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>
#include <aleph/topology/io/SimplicialComplexReader.hh>

// Persistence diagrams ------------------------------------------------
//
// Persistence diagram class, various norms, and conversion functions,
// as well as the standard distance calculations plus kernels.

#include <aleph/persistenceDiagrams/Calculation.hh>
#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/PersistenceIndicatorFunction.hh>

#include <aleph/persistenceDiagrams/distances/Bottleneck.hh>
#include <aleph/persistenceDiagrams/distances/Hausdorff.hh>
#include <aleph/persistenceDiagrams/distances/Wasserstein.hh>

#include <aleph/persistenceDiagrams/kernels/MultiScaleKernel.hh>

#include <aleph/persistenceDiagrams/io/Raw.hh>

// Persistent homology calculation algorithms --------------------------
//
// The "general" calculation is being exposed as well as specific
// classes for connected components.

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/ConnectedComponents.hh>
#include <aleph/persistentHomology/PersistencePairing.hh>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <tuple>

#include <cmath>

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
using Representation     = aleph::defaults::Representation;
using ReductionAlgorithm = aleph::defaults::ReductionAlgorithm;

// Select a default value of calculating nearest neighbours. This is
// only relevant for those functions that create complexes from data
// that is unstructured.
template <class Distance> using NearestNeighbours
#ifdef ALEPH_WITH_FLANN
  = aleph::geometry::FLANN<PointCloud, Distance>;
#else
  = aleph::geometry::BruteForce<PointCloud, Distance>;
#endif

void wrapSimplex( py::module& m )
{
  py::class_<Simplex>(m, "Simplex")
    .def( py::init<>() )
    .def( py::init<VertexType, DataType>() )
    .def( py::init<Simplex, DataType>() )
    .def( py::init(
            [] ( py::list vertices_ )
            {
              std::vector<VertexType> vertices;
              for( auto handle : vertices_ )
                vertices.push_back( py::cast<VertexType>( handle ) );

              return new Simplex( vertices.begin(), vertices.end() );
            }
          )
    )
    .def( py::init(
            [] ( py::list vertices_, DataType dataType )
            {
              std::vector<VertexType> vertices;
              for( auto handle : vertices_ )
                vertices.push_back( py::cast<VertexType>( handle ) );

              return new Simplex( vertices.begin(), vertices.end(), dataType );
            }
          )
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
    .def( py::init(
            [] ( py::list simplices_ )
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

              return new SimplicialComplex( simplices.begin(), simplices.end() );
            }
          )
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
      [] ( PersistenceDiagram& D )
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
    )
    .def( "__array__",
      [] ( PersistencePairing& pairing )
      {
        auto n_points    = pairing.size();
        DataType* buffer = new DataType[ 2*n_points ];

        for( auto it = pairing.begin(); it != pairing.end(); ++it )
        {
          auto index        = std::distance( pairing.begin(), it );
          buffer[2*index  ] = it->first;
          buffer[2*index+1] = it->second;
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

  // Wraps a function that calculates a zero-dimensional persistence
  // diagram along with a persistence pairing. This will permit us a
  // simpler assignment of total persistence values to layers.
  m.def( "calculateZeroDimensionalPersistenceDiagramAndPairing",
    [] ( const SimplicialComplex& K, DataType unpairedData )
    {
      // Tells the calculation procedure that a pairing should be
      // calculated alongside the persistence diagram.
      using Traits
        = aleph::traits::PersistencePairingCalculation<PersistencePairing>;

      using Point = typename PersistenceDiagram::Point;
      auto tuple  = aleph::calculateZeroDimensionalPersistenceDiagram<Simplex, Traits>( K );
      auto&& pd   = std::get<0>( tuple );

      if( std::isfinite( unpairedData ) )
      {
        std::transform( pd.begin(), pd.end(), pd.begin(),
          [&unpairedData] ( const Point& p )
          {
            if( p.isUnpaired() )
            {
              return Point( p.x(), unpairedData );
            }
            else
              return Point( p );
          }
        );
      }

      return tuple;
    },
    py::arg("K"),
    py::arg("unpairedData") = std::numeric_limits<DataType>::infinity()
  );

  // Wraps a function that calculates a zero-dimensional persistence
  // diagram from a `numpy` matrix. The matrix is assumed to contain
  // edge weights of a bipartite graph. Depending on the parameters,
  // the graph will be filtered from large weights to small ones, or
  // vice versa.
  m.def( "calculateZeroDimensionalPersistenceDiagramOfMatrix",
    [] ( py::array_t<double> M, bool reverseFiltration, DataType vertexWeight, DataType unpairedData )
    {
      py::buffer_info bufferInfo = M.request();

      if( bufferInfo.ndim != 2 || bufferInfo.shape.size() != 2 )
        throw std::runtime_error( "Only two-dimensional buffers are supported" );

      if( bufferInfo.format != py::format_descriptor<DataType>::format() )
        throw std::runtime_error( "Buffer format is not consistent with data type" );

      auto n = bufferInfo.shape[0];
      auto m = bufferInfo.shape[1];

      std::vector<Simplex> simplices;
      simplices.reserve( static_cast<std::size_t>( (n+m) + (n*m) ) );

      for( VertexType v = 0; v < VertexType(n+m); v++ )
        simplices.push_back( Simplex( VertexType(v), vertexWeight) );

      // Determine the proper stride for accessing the array. While it
      // is very probable that this just defaults to the C array index
      // that is so common, i.e. i*m + j, I want to be sure.
      std::size_t rowStride = std::size_t( bufferInfo.strides[0] ) / sizeof(DataType);
      std::size_t colStride = std::size_t( bufferInfo.strides[1] ) / sizeof(DataType);

      for( VertexType u = 0; u < VertexType(n); u++ )
      {
        for( VertexType v = 0; v < VertexType(m); v++ )
        {
          simplices.push_back(
            Simplex(
              {
                u,
                VertexType(v+n)
              },
              reinterpret_cast<DataType*>( bufferInfo.ptr )[u*rowStride+v*colStride] ) );
        }
      }

      SimplicialComplex K( simplices.begin(), simplices.end() );

      if( reverseFiltration )
      {
        K.sort( aleph::topology::filtrations::Data<Simplex,
                                                   std::greater<DataType> >()
        );
      }
      else
        K.sort( aleph::topology::filtrations::Data<Simplex>() );

      using Point = typename PersistenceDiagram::Point;
      auto tuple  = aleph::calculateZeroDimensionalPersistenceDiagram<Simplex>( K );
      auto&& pd   = std::get<0>( tuple );

      if( std::isfinite( unpairedData ) )
      {
        std::transform( pd.begin(), pd.end(), pd.begin(),
          [&unpairedData] ( const Point& p )
          {
            if( p.isUnpaired() )
            {
              return Point( p.x(), unpairedData );
            }
            else
              return Point( p );
          }
        );
      }

      return pd;
    },
    py::arg("M"),
    py::arg("reverseFiltration") = true,
    py::arg("vertexWeight")      = DataType(1.0),
    py::arg("unpairedData")      = std::numeric_limits<DataType>::infinity()
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
    .def( py::init(
            [] ( const PersistenceDiagram& D )
            {
              return new StepFunction( aleph::persistenceIndicatorFunction( D ) );
            }
      )
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
      "k"_a        = 2.0,
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

/**
  Provides wrapper functions for different ways of calculating
  a Vietoris--Rips complex of a data set. These functions have
  been written for the purpose of being applied to batches. No
  speedups to be expected.
*/

void wrapVietorisRipsComplexCalculation( py::module& m )
{
  using namespace pybind11::literals;

  // Wraps a function that calculates a Vietoris--Rips complex up to
  // a pre-defined dimension, returning a persistence diagram, and a
  // corresponding pairing.
  m.def( "vietoris_rips_from_matrix",
    [] ( py::array_t<double> M, unsigned max_dimension )
    {
      py::buffer_info bufferInfo = M.request();

      if( bufferInfo.ndim != 2 || bufferInfo.shape.size() != 2 )
        throw std::runtime_error( "Only two-dimensional buffers are supported" );

      if( bufferInfo.format != py::format_descriptor<DataType>::format() )
        throw std::runtime_error( "Buffer format is not consistent with data type" );

      auto n = bufferInfo.shape[0];
      auto m = bufferInfo.shape[1];

      // We require this to be a distance matrix, even though we do not
      if( n != m )
        throw std::runtime_error( "Unable to handle rectangular matrices" );

      std::vector<Simplex> simplices;
      simplices.reserve( static_cast<std::size_t>( n + (n * (n - 1) / 2 ) ) );

      // Create vertices following the idea of a Vietoris--Rips complex
      // that handles a distance function. This is the only sane thing,
      // for we are getting a distance matrix as an input.
      for( VertexType v = 0; v < VertexType(n); v++ )
        simplices.push_back( Simplex( VertexType(v), DataType() ) );

      // Determine the proper stride for accessing the array. While it
      // is very probable that this just defaults to the C array index
      // that is so common, i.e. i*m + j, I want to be sure.
      std::size_t rowStride = std::size_t( bufferInfo.strides[0] ) / sizeof(DataType);
      std::size_t colStride = std::size_t( bufferInfo.strides[1] ) / sizeof(DataType);

      // This loop assumes that the diagonal values of the distance
      // matrix are zero. In theory, they could even contain random
      // values because we ignore them outright.
      for( VertexType u = 0; u < VertexType(n); u++ )
      {
        for( VertexType v = u + 1; v < VertexType(n); v++ )
        {
          // Add the corresponding edge, with the weight being
          // calculated from the distance matrix.
          simplices.push_back(
            Simplex(
              {
                u, v
              },
              reinterpret_cast<DataType*>( bufferInfo.ptr )[u*rowStride+v*colStride] ) );
        }
      }

      SimplicialComplex K( simplices.begin(), simplices.end() );

      // Perform the desired expansion in case it has been requested by
      // the client.
      if( max_dimension != 0 )
      {
        RipsExpander ripsExpander;
        K = ripsExpander( K, max_dimension );
        K = ripsExpander.assignMaximumWeight( K );
      }

      // Bring the complex into proper filtration order, following the
      // idea of a distance filtration.
      K.sort( aleph::topology::filtrations::Data<Simplex>() );

      auto boundaryMatrix = aleph::topology::makeBoundaryMatrix<Representation>( K );
      auto pairing        = aleph::calculatePersistencePairing<ReductionAlgorithm>( boundaryMatrix.dualize() );
      auto diagrams       = aleph::makePersistenceDiagrams( pairing, K );

      return diagrams;
    },
    py::arg("M"),
    py::arg("max_dimension") = 0
  );

  // Wraps a function that calculates a Vietoris--Rips complex up to
  // dimension 1 (not including cycles), and returns a corresponding
  // persistence pairing. The pairing will consist of the destroyer,
  // i.e. the edge that merges a connected component.
  m.def( "vietoris_rips_from_matrix_1d",
    [] ( py::array_t<double> M )
    {
      py::buffer_info bufferInfo = M.request();

      if( bufferInfo.ndim != 2 || bufferInfo.shape.size() != 2 )
        throw std::runtime_error( "Only two-dimensional buffers are supported" );

      if( bufferInfo.format != py::format_descriptor<DataType>::format() )
        throw std::runtime_error( "Buffer format is not consistent with data type" );

      auto n = bufferInfo.shape[0];
      auto m = bufferInfo.shape[1];

      // We require this to be a distance matrix, even though we do not
      if( n != m )
        throw std::runtime_error( "Unable to handle rectangular matrices" );

      std::vector<Simplex> simplices;
      simplices.reserve( static_cast<std::size_t>( n + (n * (n - 1) / 2 ) ) );

      // Create vertices following the idea of a Vietoris--Rips complex
      // that handles a distance function. This is the only sane thing,
      // for we are getting a distance matrix as an input.
      for( VertexType v = 0; v < VertexType(n); v++ )
        simplices.push_back( Simplex( VertexType(v), DataType() ) );

      // Determine the proper stride for accessing the array. While it
      // is very probable that this just defaults to the C array index
      // that is so common, i.e. i*m + j, I want to be sure.
      std::size_t rowStride = std::size_t( bufferInfo.strides[0] ) / sizeof(DataType);
      std::size_t colStride = std::size_t( bufferInfo.strides[1] ) / sizeof(DataType);

      // This loop assumes that the diagonal values of the distance
      // matrix are zero. In theory, they could even contain random
      // values because we ignore them outright.
      for( VertexType u = 0; u < VertexType(n); u++ )
      {
        for( VertexType v = u + 1; v < VertexType(n); v++ )
        {
          // Add the corresponding edge, with the weight being
          // calculated from the distance matrix.
          simplices.push_back(
            Simplex(
              {
                u, v
              },
              reinterpret_cast<DataType*>( bufferInfo.ptr )[u*rowStride+v*colStride] ) );
        }
      }

      // Bring the complex into proper filtration order, following the
      // idea of a distance filtration.
      SimplicialComplex K( simplices.begin(), simplices.end() );
      K.sort( aleph::topology::filtrations::Data<Simplex>() );

      auto boundaryMatrix = aleph::topology::makeBoundaryMatrix<Representation>( K );
      auto pairing        = aleph::calculatePersistencePairing<ReductionAlgorithm>( boundaryMatrix.dualize() );

      using Pair = decltype( pairing )::ValueType;

      // This follows the terminology of a down-stream task, in which we
      // think of this as a (very convoluted!) process of selecting edge
      // subsets, i.e. those that belong to the MST.
      std::vector<Pair> selected_edges;
      selected_edges.reserve( static_cast<std::size_t>( n ) );

      for( auto&& pair : pairing )
      {
        // Do not add a pair for the unpaired vertex, as it will result
        // in a cycle edge.
        if( pair.second > K.size() )
          continue;

        // Extract the proper vertex indices of the edge in order to
        // make this a proper selection of edges.
        else
        {
          auto edge = K[ pair.second ];
          selected_edges.push_back( std::make_pair( edge[1], edge[0] ) );
        }
      }

      return selected_edges;
    },
    py::arg("M")
  );

  // Wraps a function that calculates a Vietoris--Rips complex up to
  // dimension 2 (including cycles), returning the persistence pairs
  // in a special format:
  //
  // A tuple of two lists will be returned. The first contains every
  // selected edge as in the 1D case. The second one contains a list
  // of 4 indices [i, j, k, l]. The first two indices describe edges
  // that *create* a cycle, while the last two indices contain *the*
  // unique (!) edge of the triangle that destroys the cycle.
  m.def( "vietoris_rips_from_matrix_2d",
    [] ( py::array_t<double> M )
    {
      py::buffer_info bufferInfo = M.request();

      if( bufferInfo.ndim != 2 || bufferInfo.shape.size() != 2 )
        throw std::runtime_error( "Only two-dimensional buffers are supported" );

      if( bufferInfo.format != py::format_descriptor<DataType>::format() )
        throw std::runtime_error( "Buffer format is not consistent with data type" );

      auto n = bufferInfo.shape[0];
      auto m = bufferInfo.shape[1];

      // We require this to be a distance matrix, even though we do not
      if( n != m )
        throw std::runtime_error( "Unable to handle rectangular matrices" );

      std::vector<Simplex> simplices;
      simplices.reserve( static_cast<std::size_t>( n + (n * (n - 1) / 2 ) ) );

      // Create vertices following the idea of a Vietoris--Rips complex
      // that handles a distance function. This is the only sane thing,
      // for we are getting a distance matrix as an input.
      for( VertexType v = 0; v < VertexType(n); v++ )
        simplices.push_back( Simplex( VertexType(v), DataType() ) );

      // Determine the proper stride for accessing the array. While it
      // is very probable that this just defaults to the C array index
      // that is so common, i.e. i*m + j, I want to be sure.
      std::size_t rowStride = std::size_t( bufferInfo.strides[0] ) / sizeof(DataType);
      std::size_t colStride = std::size_t( bufferInfo.strides[1] ) / sizeof(DataType);

      // This loop assumes that the diagonal values of the distance
      // matrix are zero. In theory, they could even contain random
      // values because we ignore them outright.
      for( VertexType u = 0; u < VertexType(n); u++ )
      {
        for( VertexType v = u + 1; v < VertexType(n); v++ )
        {
          // Add the corresponding edge, with the weight being
          // calculated from the distance matrix.
          simplices.push_back(
            Simplex(
              {
                u, v
              },
              reinterpret_cast<DataType*>( bufferInfo.ptr )[u*rowStride+v*colStride] ) );
        }
      }

      SimplicialComplex K( simplices.begin(), simplices.end() );

      // Perform the Vietoris--Rips expansion to be able to handle
      // cycles in addition to 1D topology.
      {
        RipsExpander ripsExpander;
        K = ripsExpander( K, 2 );
        K = ripsExpander.assignMaximumWeight( K );
      }

      // Bring the complex into proper filtration order, following the
      // idea of a distance filtration.
      K.sort( aleph::topology::filtrations::Data<Simplex>() );

      auto boundaryMatrix = aleph::topology::makeBoundaryMatrix<Representation>( K );
      auto pairing        = aleph::calculatePersistencePairing<ReductionAlgorithm>( boundaryMatrix.dualize() );

      using Index = decltype( pairing )::IndexType;
      using Pair  = decltype( pairing )::ValueType;

      // This follows the terminology of a down-stream task, in which we
      // think of this as a (very convoluted!) process of selecting edge
      // subsets, i.e. those that belong to the MST.
      std::vector<Pair> selected_edges;
      selected_edges.reserve( static_cast<std::size_t>( n ) );

      // Describes a pair of edges: the first two indices describe the
      // indices of the creator edge, while the last two correspond to
      // the indices of the destroyer edge.
      using EdgePair = std::tuple<Index, Index, Index, Index>;

      std::vector<EdgePair> selected_pairs;
      selected_pairs.reserve( static_cast<std::size_t>( n * (n - 1) / 2 ) );

      for( auto&& pair : pairing )
      {
        auto&& creator = K[ pair.first ];

        // Handle 0-dimensional features as before
        if( creator.dimension() == 0 )
        {
          // Do not add a pair for the unpaired vertex, as it will result
          // in a cycle edge.
          if( pair.second > K.size() )
            continue;

          // Extract the proper vertex indices of the edge in order to
          // make this a proper selection of edges.
          else
          {
            auto u = K[ pair.first ][0]; // creator is always a vertex
            auto e = K[ pair.second ];   // destroyer is always an edge

            // Select the proper vertex index for the destroyer; this is
            // just the *other* vertex of the edge.
            auto v = e[0] != u ? e[0] : e[1];

            selected_edges.push_back( std::make_pair( u, v ) );
          }
        }

        // Handle 1-dimensional features by searching for the maximum
        // edge in the destroying triangle.
        else if( creator.dimension() == 1 )
        {
          auto&& destroyer = K[ pair.second ];

          // Triangle (u, v, w) gives rise to the following edges:
          // - [u, v] (d1)
          // - [u, w] (d2)
          // - [v, w] (d3)
          //
          // The way we access them ensures some lexicographical
          // ordering already, so we do not have to sort them in
          // a later stage.
          auto u  = destroyer[2];
          auto v  = destroyer[1];
          auto w  = destroyer[0];

          auto d1 = reinterpret_cast<DataType*>( bufferInfo.ptr )[u*rowStride+v*colStride];
          auto d2 = reinterpret_cast<DataType*>( bufferInfo.ptr )[u*rowStride+w*colStride];
          auto d3 = reinterpret_cast<DataType*>( bufferInfo.ptr )[v*rowStride+w*colStride];

          // Prepare creator; there's not much we can do here except for
          // returning the edges in the right order (hence swapped index
          // access).

          auto i = creator[1];
          auto j = creator[0];

          // For the destroyer, we have to pick the edge with the
          // largest weight, according to {d1, d2, d3}. We do not
          // know which edge is going to be the largest, so we do
          // a check for all of them.
          //
          // Notice that the ordering of the if clauses implies a
          // sane ordering in case all of the edges have the same
          // weight.

          auto k = i;
          auto l = j;

          if( d1 >= creator.data() )
          {
            k = u; 
            l = v;
          }
          else if( d2 >= creator.data() )
          {
            k = u;
            l = w;
          }
          else if( d3 >= creator.data() )
          {
            k = v;
            l = w;
          }

          selected_pairs.push_back( std::make_tuple(i, j, k, l) );
        }
      }

      return std::make_pair( selected_edges, selected_pairs );
    },
    py::arg("M")
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
  wrapVietorisRipsComplexCalculation(m);
}
