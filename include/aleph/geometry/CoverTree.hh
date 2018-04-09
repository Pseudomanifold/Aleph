#ifndef ALEPH_GEOMETRY_COVER_TREE_HH__
#define ALEPH_GEOMETRY_COVER_TREE_HH__

#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

// FIXME: remove after debugging
#include <iostream>

namespace aleph
{

namespace geometry
{

/**
  @class CoverTree
  @brief Generic cover tree data structure

  This class models a cover tree data structure, as described in the
  paper "Cover trees for nearest neighbor" by Beygelzimer et al.; an
  implementation is given at:

    http://hunch.net/~jl/projects/cover_tree/cover_tree.html

  This implementation attempts to be as generic as possible.
*/

template <class Point, class Metric> class CoverTree
{
public:

  class Node
  {
  public:

    /** Creates a new node that stores a point */
    Node( const Point& point )
      : _point( point )
    {
    }

    /** The point stored in the node */
    Point _point;

    /**
      All children of the node. Their order depends on the insertion
      order into the data set.
    */

    std::vector< std::unique_ptr<Node> > _children;
  };

  /**
    Creates a cover tree from a range of points. The client can specify
    whether to use batch insertion or standard insertion.

    @param begin Input iterator to begin of point range
    @param end   Input iterator to end of point range
    @param batch Flag indicating whether batch insertion should be performed
  */

  template <class InputIterator> CoverTree(
    InputIterator begin, InputIterator end,
    bool batch = false )
  {
    if( batch )
      throw std::runtime_error( "Not yet implemented" );

    auto n = std::distance( begin, end );
    if( n == 0 )
      throw std::runtime_error( "Cover tree must be non-empty" );

    _root = std::unique_ptr<Node>( new Node( *begin ) );

    for( auto it = std::next( begin ); it != end; ++it )
    {
      this->insert( *it,
                    { _root.get() },
                    std::numeric_limits<unsigned>::max() );
    }
  }

private:

  /**
    Auxiliary function for calculating the distance between a point and
    a given set of nodes.
  */

  double distanceToSet( const Point& p, const std::vector<Node*>& Q )
  {
    Metric metric;

    double d = std::numeric_limits<double>::max();
    for( auto&& q : Q )
      d = std::min( d, double( metric(p, q->_point) ) );

    return d;
  }

  bool insert( const Point& p, const std::vector<Node*>& Qi, unsigned i )
  {
    std::cerr << "Inserting point " << p << " at level " << i << "\n";

    Metric metric;

    // Contains all nodes that satisfy the distance criterion for this
    // level.
    std::vector<Node*> Qj;

    double d = std::numeric_limits<double>::max();
    for( auto&& Q : Qi )
    {
      for( auto&& child : Q->_children )
      {
        double distance = double( metric( p, child->_point ) );
        if( unsigned( std::log2( distance ) ) <= i )
        {
          d = std::min( d, distance );
          Qj.push_back( child.get() );
        }
      }
    }

    d = std::log2( d );

    if( d > i )
      return false;
    else
    {
      auto parentFound = insert( p, Qj, i-1 );
      auto distance    = unsigned( std::log2( distanceToSet(p, Qi) ) );

      if( !parentFound && distance < i )
        return true;
      else
        return false;
    }
  }

  /** Root pointer of the tree */
  std::unique_ptr<Node> _root;
};

} // namespace geometry

} // namespace aleph

#endif
