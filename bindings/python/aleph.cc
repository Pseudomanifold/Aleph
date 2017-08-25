#include <pybind11/pybind11.h>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

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
    .def( "__repr__",
          [] ( const Simplex& simplex )
          {
            std::ostringstream stream;
            stream << simplex;

            return stream.str();
          }
    );
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
                auto&& vertices_ = py::cast<py::list>( simplexHandle );

                std::vector<VertexType> vertices;
                for( auto vertexHandle : vertices_ )
                  vertices.push_back( py::cast<VertexType>( vertexHandle ) );

                simplices.push_back( Simplex( vertices.begin(), vertices.end() ) );
              }

            }

            new (&instance) SimplicialComplex( simplices.begin(), simplices.end() );
          }
    )
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
    );
}

PYBIND11_PLUGIN(aleph)
{
  py::module m("aleph", "Python bindings for Aleph, a library for exploring persistent homology");

  wrapSimplex(m);
  wrapSimplicialComplex(m);

  return m.ptr();
}
