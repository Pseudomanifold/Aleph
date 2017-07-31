#ifndef ALEPH_PERSISTENCE_DIAGRAMS_IO_JSON_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_IO_JSON_HH__

#include <fstream>
#include <ostream>
#include <string>

namespace aleph
{

namespace io
{

/**
  Writes a persistence diagram to an output stream, using the JSON
  format. The diagram will be serialized such that its points will
  be stored in the field ``data`` as a two-dimensional array. Note
  that infinite values will be encoded as strings. Additional data
  about the diagram, e.g. its dimension, are stored in name--value
  pairs.
*/

template <class Diagram> void writeJSON( std::ostream& o, const Diagram& D, const std::string& name = std::string() )
{
  std::string level = "  ";

  o << "{\n";

  o << level << "\"betti\": "     << D.betti()     << ",\n"
    << level << "\"dimension\": " << D.dimension() << ",\n";

  if( !name.empty() )
    o << level << "\"name\": " << "\"" << name << "\",\n";

  o << level << "\"size\": "      << D.size()      << ",\n"
    << level << "\"diagram\": "   << "[\n";

  for( auto it = D.begin(); it != D.end(); ++it )
  {
    if( it != D.begin() )
      o << ",\n";

    o << level << level << "["
                        << "\"" << it->x() << "\""
                        << ","
                        << "\"" << it->y() << "\""
                        << "]";
  }

  o << "\n"
    << level << "]\n"
    << "}";
}

} // namespace io

} // namespace aleph

#endif
