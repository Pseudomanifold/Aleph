#include <pybind11/pybind11.h>

#include <aleph/topology/Simplex.hh>

namespace py = pybind11;

using DataType   = double;
using VertexType = unsigned;

using Simplex    = aleph::topology::Simplex<DataType, VertexType>;

PYBIND11_PLUGIN(aleph)
{
  py::module m("aleph", "Python bindings for Aleph, a library for exploring persistent homology");

  py::class_<Simplex>(m, "Simplex")
    .def( py::init<>() )
    .def( py::init<VertexType, DataType>() )
    .def( py::init<Simplex, DataType>() )
    .def( "__init__",
          [] ( Simplex& instance, py::list vertices_, DataType dataType )
          {
            std::vector<VertexType> vertices;
            for( auto vertexHandle : vertices_ )
              vertices.push_back( py::cast<VertexType>( vertexHandle ) );

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

  return m.ptr();
}
