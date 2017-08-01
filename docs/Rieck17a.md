---
title: Hierarchies and Ranks for Persistence Pairs
---

# Hierarchies and Ranks for Persistence Pairs

In this paper, which is available as a [preprint on my personal
homepage](http://bastian.rieck.ru/research/TopoInVis2017_Hierarchies.pdf),
we developed a novel hierarchy, the *interlevel set persistence
hierarchy* that relates the critical pairs stored in a persistence
diagram with each other.

Aleph ships with a tool called `interlevel_set_persistence_hierarchy`
that permits the calculation of this hierarchy for 1D data and 2D VTK
data. The tool supports the calculation of both superlevel and sublevel
sets as explained in the paper.

## Simple 1D example

As a first simple 1D example, let us compare the two functions described
in the paper. We first need the input data. Put the following two lines
in a file `1D.txt`:

    3 1 6 5 8 2 7 4
    3 1 8 2 7 5 6 4

Each line represents a function, or rather, its *image*. Notice that the
domain of the function does not matter for topological calculations.

New, we call the tool and take a look at the hierarchy based on sublevel
sets:


    $ ./interlevel_set_persistence_hierarchy 1D.txt

Among others, the output should contain the following:

    0: 1    inf
    1: 5    6
    2: 2    8
    3: 4    7

    0 -- 1
    2 -- 3
    0 -- 2


    0: 1    inf
    1: 2    8
    2: 5    6
    3: 4    7

    3 -- 2
    1 -- 3
    0 -- 1

The mapping shows the index that was assigned to a given critical point.
The remaining lines indicate how those points are connected. Individual
hierarchies are separated by two newlines, so they can be parsed easily
later on.

Putting the first hierarchy into a file `01.txt` and the second
hierarchy into a file `02.txt`, we can use another script to compare
them:

    $ ./compare_interlevel_set_persistence_hierarchies.py 01.txt 02.txt

The output of this script is unfortunately still somewhat rough&nbsp;(if
you want to help improve the script, please consider [contributing to
Aleph](https://github.com/Submanifold/Aleph/blob/master/CONTRIBUTING.md)).
The file `/tmp/M.txt` should now contain the following output:

    0.000000000000000000e+00 2.000000000000000000e+03
    2.000000000000000000e+03 0.000000000000000000e+00

The off-diagonal values indicate that the two hierarchies are different.
