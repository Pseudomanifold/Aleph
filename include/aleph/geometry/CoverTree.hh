#ifndef ALEPH_GEOMETRY_COVER_TREE_HH__
#define ALEPH_GEOMETRY_COVER_TREE_HH__

#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <queue>
#include <stack>
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
    Node( const Point& point, long level )
      : _point( point )
      , _level( level )
    {
    }

    /** Calculates current covering distance of the node */
    double coveringDistance() const noexcept
    {
      return std::pow( coveringConstant, static_cast<double>( _level ) );
    }

    /** Calculates current separating distance of the node */
    double separatingDistance() const noexcept
    {
      return std::pow( coveringConstant, static_cast<double>( _level - 1 ) );
    }

    /** @returns true if the node is a leaf node */
    bool isLeaf() const noexcept
    {
      return _children.empty();
    }

    void insert( const Point& p )
    {
      auto d = Metric()( _point, p );

      std::cerr << __FUNCTION__ << ": Inserting " << p << "\n";
      std::cerr << __FUNCTION__ << ": Covering distance           = " << this->coveringDistance() << "\n";
      std::cerr << __FUNCTION__ << ": Distance from point to root = " << d << "\n";

      if( d > this->coveringDistance() )
      {
        while( d > 2 * this->coveringDistance() )
        {
          std::cerr << __FUNCTION__ << ": Distance is bigger than covering distance; need to raise level of tree\n";

          // -----------------------------------------------------------
          //
          // Find a leaf node that can become the new root node with
          // a raised level. If the tree only contains the root node
          // its level can be adjusted multiple times.

          std::stack<Node*> nodes;
          nodes.push( this );

          Node* leaf   = nullptr;
          Node* parent = nullptr;

          // Special case: the root itself is a leaf node; this happens
          // at the beginning of the insertion process and means that a
          // level adjustment has to be performed.
          if( this->isLeaf() )
          {
            this->_level += 1;
            continue;
          }

          while( !nodes.empty() )
          {
            auto node = nodes.top();
            nodes.pop();

            for( auto&& child : node->_children )
            {
              if( child->isLeaf() )
              {
                leaf   = child.get();
                parent = node;
                break;
              }
              else
                nodes.push( child.get() );
            }
          }

          std::cerr << __FUNCTION__ << ": Found leaf node " << leaf->_point << "\n";

          // There is no leaf, so there is nothing to do and we just
          // skip to the bottom where we add the current node as the
          // new root of the tree.
          if( !leaf )
          {
            std::cerr << __FUNCTION__ << ": Unable to identify leaf node\n";
            break;
          }

          assert( leaf );
          assert( parent );

          // Remove leaf from subtree ----------------------------------
          //
          // The previous tree does not contain the leaf node any more,
          // and it can be added as the new root node.

          std::unique_ptr<Node> leaf_ptr( nullptr );

          parent->_children.erase(
            std::remove_if(
              parent->_children.begin(), parent->_children.end(),
              [&leaf, &leaf_ptr] ( std::unique_ptr<Node>& child )
              {
                if( child.get() == leaf )
                {
                  leaf_ptr = std::move( child );
                  return true;
                }

                return false;
              }
            ),
            parent->_children.end()
          );

          // Make leaf node the new root node --------------------------
          //
          // Notice that this does increase the level of the current
          // root. This has to be adjusted for.

          auto oldRoot
            = std::unique_ptr<Node>( new Node( this->_point, this->_level ) );

          for( auto&& child : _children )
            oldRoot->_children.push_back( std::move( child ) );

          _point = leaf->_point;
          _level = _level + 1;

          _children.clear();
          _children.push_back( std::move( oldRoot ) );

          // Since the root of the tree changed, we also have to update
          // the distance calculation.
          d = Metric()( _point, p );

          std::cerr << __FUNCTION__ << ": " << leaf->_point << " is the new root\n";
        }

        // Make current point the new root -----------------------------
        //
        // So far, the new point has not yet been inserted into the
        // tree. This needs to be done now while the cover is valid
        // again.

        auto oldRoot
          = std::unique_ptr<Node>( new Node( this->_point, this->_level ) );

        for( auto&& child : _children )
          oldRoot->_children.push_back( std::move( child ) );

        _point = p;
        _level = _level + 1;

        _children.clear();
        _children.push_back( std::move( oldRoot ) );

        return;
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

      // Add the new point as a child of the current root node. Note the
      // level adjustment.
      _children.push_back( std::unique_ptr<Node>( new Node( p, _level - 1 ) ) );
    }

    Point _point; //< The point stored in the node
    long  _level; //< The level of the node

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
      _root = std::unique_ptr<Node>( new Node(p,0) );
    else
      _root->insert( p );
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

  // Validity checks ---------------------------------------------------
  //
  // These are called by debug code and tests to ensure that the cover
  // tree is correct.

  /**
    Checks the level invariant in the cover tree. The level invariant
    states that the level \f$l$ of the *direct* child of a node \f$p$
    is the level \f$l'-1$, where \f$l'$ is the level of \f$p$.
  */

  bool checkLevelInvariant() const noexcept
  {
    std::queue<const Node*> nodes;
    nodes.push( _root.get() );

    // Initialize the level artificially in order to simplify the code
    // below.
    long level = _root->_level + 1;

    while( !nodes.empty() )
    {
      auto n = nodes.size();

      for( decltype(n) i = 0; i < n; i++ )
      {
        auto&& node = nodes.front();

        if( node->_level != level - 1 )
          return false;

        for( auto&& child : node->_children )
          nodes.push( child.get() );

        nodes.pop();
      }

      level -= 1;
    }

    return true;
  }

  /**
    Checks the covering invariant in the cover tree. The covering
    invariant states that the distance between a child and parent
    node is bounded by the covering distance.
  */

  bool checkCoveringInvariant() const noexcept
  {
    std::queue<const Node*> nodes;
    nodes.push( _root.get() );

    while( !nodes.empty() )
    {
      auto n = nodes.size();
      for( decltype(n) i = 0; i < n; i++ )
      {
        auto&& parent = nodes.front();

        for( auto&& child : parent->_children )
        {
          auto d = Metric()( parent->_point, child->_point );
          if( d > parent->coveringDistance() )
          {
            std::cerr << __FUNCTION__ << ": Covering invariant is violated by ("
                      << parent->_point << "," << child->_point << "): "
                      << d << " > " << parent->coveringDistance() << "\n";

            return false;
          }

          nodes.push( child.get() );
        }

        // All children of the current parent node have been processed,
        // so we can remove it.
        nodes.pop();
      }
    }

    return true;
  }

  /**
    Checks the separating invariant in the cover tree. The separating
    invariant states that the distance between a child and its parent
    is larger than the separating distance.
  */

  bool checkSeparatingInvariant() const noexcept
  {
    std::queue<const Node*> nodes;
    nodes.push( _root.get() );

    while( !nodes.empty() )
    {
      auto n = nodes.size();
      for( decltype(n) i = 0; i < n; i++ )
      {
        auto&& parent = nodes.front();

        for( auto it1 = parent->_children.begin(); it1 != parent->_children.end(); ++it1 )
        {
          for( auto it2 = std::next(it1); it2 != parent->_children.end(); ++it2 )
          {
            // The distance between the two points must by necessity be
            // larger than the separating distance of their parent node
            // in order to satisfy the separating property.

            auto&& p = (*it1)->_point;
            auto&& q = (*it2)->_point;
            auto d   = Metric()(p, q);

            if( d <= parent->separatingDistance() )
            {
              std::cerr << __FUNCTION__ << ": Separating invariant is violated by ("
                        << p << "," << q << "): "
                        << d << " > " << parent->separatingDistance() << "\n";

              return false;
            }
          }

          // Add the child such that the next level of the tree can be
          // visited.
          nodes.push( it1->get() );
        }

        // All children of the current parent node have been processed,
        // so we can remove it.
        nodes.pop();
      }
    }

    return true;
  }

  /**
    Generic validity check of the tree. Combines *all* validity
    criteria, i.e. the level invariant, the covering invariant,
    and the separating invariant.
  */

  bool isValid() const noexcept
  {
    return    this->checkLevelInvariant()
           && this->checkCoveringInvariant()
           && this->checkSeparatingInvariant();
  }

private:

  /** Root pointer of the tree */
  std::unique_ptr<Node> _root;
};

} // namespace geometry

} // namespace aleph

#endif
