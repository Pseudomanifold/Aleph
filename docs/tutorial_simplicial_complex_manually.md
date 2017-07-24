---
title: Building a simplicial complex manually
---

# Building a simplicial complex manually

In contrast to other libraries for topological data analysis, Aleph is
built around a set of classes that facilitate the description of
topological objects.  The most relevant classes in this context are
`Simplex` and `SimplicialComplex`.  In practice, you probably never have
to build your own simplicial complex fully from scratch&mdash;we build
one here in order to illustrate the concepts.

## Building a simplex

The `Simplex` class models a simplex of arbitrary dimensions that
contains a *weight*, or&mdash;more generically&mdash;a *data* value, and
a list of *vertices*. The class is templated to support different choice
for both the data type and the vertex type. Let us start by building
some simplices:

```cpp
using namespace aleph::topology;

Simplex s1( {1}, 1 );
Simplex s2( {2}, 2 );
Simplex s21( {2,1}, 3 );

// Output: 1,2,3
std::cout << s1.data() << "," << s2.data() << "," << s21.data() << "\n";

// Data can also be changed
s21.setData( 42 );
```

The `Simplex` class offers some iterators as well. Vertices will be
traversed in reverse order. Faces will be traversed in lexicographical
order:

```cpp
// Output:
// - 2
// - 1
for( auto&& vertex : s21 )
  std::cout << "- " << vertex << "\n";

// Output:
// - {1}
// - {2}
for( auto it = s21.begin_boundary(); it != s21.end_boundary(); ++it )
  std::cout << "- " << *it << "\n";
```

As you can see from this example, the `Simplex` class also supports
output operations. Simplices will be denoted by curly braces, with the
data value only being shown if it is non-zero. Note that the boundary
iterator is incapable of displaying weights because only the weight for
the co-face, i.e. the original simplex, has been specified.

There are all sorts of other queries for simplices:

```cpp
std::cout << s1.dimension()  << "\n"; // 0
std::cout << s2.dimension()  << "\n"; // 0
std::cout << s21.dimension() << "\n"; // 1
std::cout << s21.size()      << "\n"; // 2 (number of vertices)
std::cout << s21.empty()     << "\n"; // false
std::cout << s21 == s1       << "\n"; // false
std::cout << s1 < s21        << "\n"; // true (lexicographical ordering)
std::cout << s21[0]          << "\n"; // 2
std::cout << s21[1]          << "\n"; // 1
std::cout << s21[2]          << "\n"; // exception
```

Simplices are also convertible to `bool`, which facilitates their use in
some algorithms. For example, the following is a valid query:

```cpp
void f( const Simplex& s )
{
  if( s )
    std::cout << "This simplex is non-empty!\n";
}
```

Please refer to the [source
code](https://github.com/Submanifold/Aleph/blob/master/include/aleph/topology/Simplex.hh)
or the API documentation for more information.

## Building a simplicial complex

The `SimplicialComplex` class models the concept of a simplicial
complex. It consists of a range of simplices that may be accessed using
different *views*. Most of the functionality is provided for
convenience purposes only&mdash;in theory, it would often be sufficient
to just provide an `std::vector` of simplices. However, the
`SimplicialComplex` class is easier to use and provides faster access to
individual simplices. Let us first build a simplicial complex that
models a triangle:

```cpp
using namespace aleph::topology;

std::vector<Simplex> simplices
  = { {1}, {2}, {4}, {2,4}, {1,2}, {1,4} };

SimplicialComplex K( simplices.begin(), simplices.end() );
```

Note that the vertex indices are not contiguous. `K` does not use vertex
indices directly, which makes it possible to store simplicial
sub-complexes and arbitrary subsets of complexes. The downside of this
approach is that the class cannot guarantee that it actually represents
a valid simplicial complexes. It is possible that certain faces are
missing from the complex. Duplicates, on the other hand, will not be
allowed automatically.

We may now iterate over `K` in numerous ways:

```cpp
// Output follows the insertion order from above. This iterator is the
// default iterator for representing a filtration.
for( auto&& simplex : K )
  std::cout << simplex << "\n";

// Output follows increasing dimensionality of simplices. This does not
// deviate from the filtration order in this example.
for( auto it = K.begin_dimension(); it != K.end_dimension(); ++it )
  std::cout << *it << "\n";

// Output follows lexicographical ordering:
// {1}
// {2}, {2,1}
// {4}, {4,1}, {4,2}
for( auto it = K.begin_lexicographical(); it != K.end_lexicographical(); ++it )
  std::cout << *it << "\n";
```

Note that these iterators are *read-only* by necessity. If you want to
replace a simplex, use the `replace()` function instead:

```cpp
for( auto it = K.begin(); it != K.end(); ++it )
{
  auto s = *it;
  s.setData( 42 );

  if( !K.replace( it, s ) )
  {
    // Simplex could not be replaced for some reason...
  }
}

```

Again, there are numerous query functions:

```cpp
std::cout << K.contains( {1,2} ) << "\n"; // true
std::cout << K.contains( {0,2} ) << "\n"; // false
std::cout << K.index( {1} )      << "\n"; // 0, because it is the first inserted simplex
std::cout << K.at( 3 )           << "\n"; // {2,4}
std::cout << K.at( 42 )          << "\n"; // exception
std::cout << K.size()            << "\n"; // 5
std::cout << K.empty()           << "\n"; // false
std::cout << K.dimension()       << "\n"; // 1, because we have at most 1-dimensional simplices
```

Finally, simplicial complexes can also be sorted. Here, we only take
a look at basic inline sorting predicates. When calculating persistent
homology, we will use predicates that are more involved:

```cpp
K.sort( [] ( const Simplex& s, const Simplex& t )
        {
          return s < t;
        }
);
```

This sorts the complex according to the lexicographical ordering induced
by the simplices.

## Loading simplicial complexes from files

The `SimplicialComplex` class is central to many calculations and
operations of Aleph. Many real-world objects, such as meshes or graphs,
can be represented as a simplicial complex. To this end, Aleph offers
several I/O classes.

We cover them in [another tutorial](tutorial_simplicial_complexes_io.md).

# More information

More information about the `SimplicialComplex` class is available in the
[source
code](https://github.com/Submanifold/Aleph/blob/master/include/aleph/topology/SimplicialComplex.hh)
or the API documentation of Aleph. We will also use simplicial
complexes when [calculating persistent homology](tutorial_persistent_homology_simplicial_complex.md).
