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

```c_cpp
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

```c_cpp

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
data value only being shown if it is non-zero.
