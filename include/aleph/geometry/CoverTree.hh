#ifndef ALEPH_GEOMETRY_COVER_TREE_HH__
#define ALEPH_GEOMETRY_COVER_TREE_HH__

#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <queue>
#include <stdexcept>
#include <vector>

// FIXME: remove after debugging
#include <iostream>

#include <cassert>
#include <cmath>

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

  This implementation attempts to be as generic as possible. It uses
  the simplified description of the cover tree, as given by Izbicki,
  Shelton in "Faster Cover Trees".
*/

template <class Point, class Metric> class CoverTree
{
public:

  /**
    Covering constant of the cover tree. It might make sense to change
    this later on in order to improve performance. Some papers set the
    constant to `1.3`.
  */

  constexpr static const double coveringConstant = 2.0;

  class Node
  {
  public:

    /** Creates a new node that stores a point */
    Node( const Point& point, unsigned level )
      : _point( point )
      , _level( level )
    {
      assert( _level >= 1 );
    }

    /** Calculates current covering distance of the node */
    double coveringDistance() const noexcept
    {
      return std::pow( coveringConstant, static_cast<double>( _level ) );
    }

    void addChild( const Point& p )
    {
      _children.push_back( std::unique_ptr<Node>( new Node( p, _level - 1 ) ) );
    }

    /** @returns true if the node is a leaf node */
    bool isLeaf() const noexcept
    {
      return _children.empty();
    }

    void insert( const Point& p )
    {
      auto d = Metric()( _point, p );

      std::cerr << __PRETTY_FUNCTION__ << ": Distance from point to root = " << d << "\n";

      if( d > this->coveringDistance() )
      {
        std::cerr << __PRETTY_FUNCTION__ << ": Distance is bigger than covering distance; need to raise level of tree\n";

        throw std::runtime_error( "Not yet implemented" );
      }

      return insert_( p );
    }

    /**
      Auxiliary function for performing the recursive insertion of a new
      node into the tree.
    */

    void insert_( const Point& p )
    {
      for( auto&& child : _children )
      {
        auto d = Metric()( child->_point, p );
        if( d <= child->coveringDistance() )
        {
          std::cerr << __PRETTY_FUNCTION__ << ": Recursive enumeration of the tree\n";

          // We found a node in which the new point can be inserted
          // *without* violating the covering invariant.
          child->insert_( p );
          return;
        }
      }

      // Add the new point as a child of the current root node. This
      // might require updating levels.

      if( _children.empty() )
        _level += 1;

      this->addChild( p );
    }

    Point    _point; //< The point stored in the node
    unsigned _level; //< The level of the node (>= 1)

    /**
      All children of the node. Their order depends on the insertion
      order into the data set.
    */

    std::vector< std::unique_ptr<Node> > _children;
  };

  /**
    Inserts a new point into the cover tree. If the tree is empty,
    the new point will become the root of the tree. Else, it shall
    be inserted according to the covering invariant.
  */

  void insert( const Point& p )
  {
    if( !_root )
      _root = std::unique_ptr<Node>( new Node(p,1) );
    else
      _root->insert( p );

    this->updateLevels();
  }

  // Pretty-printing function for the tree; this is only meant for
  // debugging purposes and could conceivably be implemented using
  // `std::ostream`.
  void print( std::ostream& o )
  {
    std::queue<const Node*> nodes;
    nodes.push( _root.get() );

    while( !nodes.empty() )
    {
      {
        auto n = nodes.size();

        for( decltype(n) i = 0; i < n; i++ )
        {
          auto&& node = nodes.front();

          if( i == 0 )
            o << node->_level << ": ";
          else
            o << " ";

          o << node->_point;

          for( auto&& child : node->_children )
            nodes.push( child.get() );

          nodes.pop();
        }

        o << "\n";
      }
    }
  }

private:

  /**
    Updates the levels of a cover tree. This function ensures that all
    information is being represented correctly. It is not mentioned in
    the original paper but appears to be necessary.
  */

  void updateLevels()
  {
    // Forward pass ----------------------------------------------------
    //
    // Determine depth of the tree. This is required in order to assign
    // levels correctly later on.

    unsigned depth = 0;

    std::queue<const Node*> nodes;
    nodes.push( _root.get() );

    while( !nodes.empty() )
    {
      auto n = nodes.size();

      for( decltype(n) i = 0; i < n; i++ )
      {
        auto node = nodes.front();
        nodes.pop();

        for( auto&& child : node->_children )
          nodes.push( child.get() );
      }

      ++depth;
    }

    assert( nodes.empty() );

    // Backward pass ---------------------------------------------------
    //
    // Set the level of the root node and update child levels. Note that
    // all nodes will be traversed because this function is dumb.

    _root->_level = depth;

    nodes.push( _root.get() );

    while( !nodes.empty() )
    {
      {
        auto&& parent = nodes.front();

        for( auto&& child : parent->_children )
        {
          assert( parent->_level >= 1 );

          child->_level = parent->_level - 1;
          nodes.push( child.get() );
        }
      }

      nodes.pop();
    }
  }

  /** Root pointer of the tree */
  std::unique_ptr<Node> _root;
};

} // namespace geometry

} // namespace aleph

#endif
