How to create the test data
===========================


Notes on recreating data in `1_2/`
----------------------------------

This test data is a CORSIKA simulation using the default HiSPARC CORSIKA
configuration. The parameters used are listed in the following table:

Option | Value
-------|-------
SEED   | 1
SEED   | 2
PRMPAR | 14
ERANGE | 1E5  1E5
THETAP | 0  0
PHIP   | 0  0

The result from the simulation `DAT000000` is converted to `corsika.h5` using
SAPPHiRE. Run the following command to convert the CORSIKA data format to HDF5:

    $ store_corsika_data 1_2/DAT000000 1_2/corsika.h5 --overwrite


Notes on recreating data in `3_4/`
----------------------------------

This test data is a CORSIKA simulation using mostly the default HiSPARC CORSIKA
configuration, except that the thinning option is also enabled. The simulation
uses the default HiSPARC input with the addition of values for the THIN option.
The parameters used are listed in the following table:

Option | Value
-------|-------
SEED   | 3
SEED   | 4
PRMPAR | 14
ERANGE | 1E5  1E5
THETAP | 0  0
PHIP   | 0  0
THIN   | 1E-4  1E10  0.5

The result from the simulation `DAT000000` is converted to `corsika.h5` using
SAPPHiRE. Run the following command to convert the CORSIKA data format to HDF5:

    $ store_corsika_data 3_4/DAT000000 3_4/corsika.h5 --overwrite


Notes on recreating `corsika_overview.h5`
-----------------------------------------

After creating the two simulations and converting them to HDF5 the overview can be created. Use the following command:

    $ generate_corsika_overview . corsika_overviw.h5`
