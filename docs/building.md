# Building Aleph

Aleph is meant to be used as a header-only library on top of which you
can develop your own projects based on persistent homology. However,
Aleph ships with numerous unit tests, some example programs, and tools
required for my current research. For building them, please clone the
repository to some local directory on your computer. Running the
following commands within this directory should be sufficient in most
cases:

    $ mkdir build
    $ cd build
    $ cmake ../
    $ make

It is advisable to test that Aleph works correctly on your system, so
you can run the unit test with:

    $ make test

Please submit any issues you may encounter.
