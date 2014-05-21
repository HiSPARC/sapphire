How to create the test data
===========================

Notes on (re)creating DIR-testdata.h5
-------------------------------------

Testdata was created by the test function itself.  It is thus a behavior
test.  So:

    1. Perform an Aires simulation
    2. Use store_aires_data.py to create a .h5 file containing ground
       particles data
    3. Perform a simulation (e.g. GroundParticlesSimulation), with the
       output group equal to SIMULATION_GROUP.
    4. Copy that file to this directory as DIR-testdata.h5
    5. Alternatively, just lift a simulation from an existing hdf5
       file.
    6. In the test file, look for lines containing:

       # For prerecording output, swap comments in following two lines

       and perform that operation.
    7. Run the test
    8. Copy the temporary test file and overwrite DIR-testdata.h5
    9. Revert the test to its original form.
    10. Run the test and verify that it succeeds.


Notes on recreating PE-testdata.h5
----------------------------------

This test data is simply downloaded raw data from the Public Database,
using the following script:

    import tables
    from datetime import datetime
    from sapphire.publicdb import download_data
    with tables.openFile('PE-testdata.h5', 'w') as data:
        download_data(data, '/s501', 501, datetime(2010, 9, 1),
                      datetime(2010, 9, 2), get_blobs=True)
