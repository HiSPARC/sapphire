How to create the test data
===========================


Notes on recreating `*_sim.h5`
------------------------------

These data files are created by running the `perform_simulation.py` script one
level up.

    $ python ../perform_simulation.py


Notes on recreating `corsika.h5`
--------------------------------

This is the same file as the `corsika.h5` file used in the CORSIKA tests for
the simulation with id `1_2`. Copy it from there.

    $ cp ../../corsika/test_data/1_2/corsika.h5 ./
