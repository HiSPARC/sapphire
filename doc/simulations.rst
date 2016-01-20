Simulation of the detection of EAS events
=========================================

Simulations are built up from several components. First the cluster needs to
be defined. This can be a set of real stations or imaginary detectors. Then
systematic detector and station errors are generated. Next a generator for
shower properties is used to determine the properties of each shower (energy
and direction). Then for each detector the shower particles that hit it are
selected. For each of these the actual detected signal strength and arrival
times (detector effects) are determined. Then the trigger conditions are
checked to see if the station triggered. If it did an event will be stored for
that station. This includes a GPS timestamp and the detector observables
(particle count and arrival time of first particle). Finally for each shower a
'coincidence' will be stored, even if it was not detected by any station. The
coincidence stores the parameters of the shower.


.. automodule:: sapphire.simulations
   :members:
   :undoc-members:

.. toctree::
   :hidden:

   simulations/base
   simulations/detector
   simulations/groundparticles
   simulations/ldf
   simulations/showerfront
