#include <pybind11/pybind11.h>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <stdexcept>

namespace py = pybind11;

using DataType   = double;
using VertexType = unsigned;

using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

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
    .def_property_readonly( "dimension", &SimplicialComplex::dimension );
}

PYBIND11_PLUGIN(aleph)
{
  py::module m("aleph", "Python bindings for Aleph, a library for exploring persistent homology");

  wrapSimplex(m);
  wrapSimplicialComplex(m);

  return m.ptr();
}
