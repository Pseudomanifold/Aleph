#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_PLUGIN(aleph)
{
  py::module m("aleph", "Python bindings for Aleph, a library for exploring persistent homology");

  return m.ptr();
}
