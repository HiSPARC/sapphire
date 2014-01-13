Notes on (re)creating testdata.h5

Testdata was created by the test function itself.  It is thus a behavior
test.


Creating new testdata.h5:

        1. Perform an Aires simulation
        2. Use store_aires_data.py to create a .h5 file containing ground
           particles data
        3. Copy that file to this directory as testdata.h5
        4. In test_groundparticles_acceptance.GroundParticlesSimulationAcceptanceTest.create_test_simulation_output
           uncomment the line with '/simulations' and comment the line
           with '/test_output'.
        5. Run the test, which crashes because of missing test_output.
        6. Copy the temporary test file (pay attention to the filename
           mentioned with "closing remaining open files") and overwrite
           testdata.h5.
        7. Revert the test to its original form.
        8. Run the test and verify that it succeeds.


Updating testdata.h5: only follow steps 4 - 8.
