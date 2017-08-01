---
title: Persistence Concepts for 2D Skeleton Evolution Analysis
---

# Persistence Concepts for 2D Skeleton Evolution Analysis

This paper, which is available as a [shortened preprint on my personal
webiste](http://bastian.rieck.ru/research/TopoInVis2017_Skeletons.pdf),
deals with&mdash;among other things&mdash;analysing time-varying data &
their topology.

This analysis relied on several statistics of persistence diagrams:

- Average persistence
- Infinity norm
- Total persistence

The tool `persistence_diagram_statistics` was developed for the purpose
of calculating these statistics and reporting them. We require a set of
persistence diagrams for demonstrating its use. Luckily, Aleph contains
[several persistence diagrams as a part of its test suite](https://github.com/Submanifold/Aleph/tree/master/tests/input),
which we may use. Having obtained them, we call the tool like this:

    $ ./persistence_diagram_statistics --invalid inf  Iris_dimension_*.txt 
    * Loading 'Iris_dimension_0.txt'...finished
    * Loading 'Iris_dimension_1.txt'...finished
    * Loading 'Iris_dimension_2.txt'...finished
    file,power,total_persistence,total_persistence_normalized,infinity_norm,average_persistence
    * Filtering all persistence pairs that contain 'inf'...
    'Iris_dimension_0.txt',2,12.02,0.0858571,0.556776,0.273387
    * Filtering all persistence pairs that contain 'inf'...
    'Iris_dimension_1.txt',2,0.0622947,0.00222481,0.104886,0.0379867
    * Filtering all persistence pairs that contain 'inf'...
    'Iris_dimension_2.txt',2,0.000452442,0.000226221,0.019248,0.0141505

The output is a CSV file (written to STDOUT initially) with the required
information in each column. The output format may change for versions of
the software, though, but the header helps us find the desired value for
each diagram easily.

Please refer to the forthcoming paper for more information. For more
details about the underlying theory and data, please visit the paper
[code repository on GitHub](https://github.com/Submanifold/Skeleton_Persistence).
