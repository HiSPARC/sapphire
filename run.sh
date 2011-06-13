#!/bin/bash -l
cd $PBS_O_WORKDIR
workon core
python detector_sim.py
