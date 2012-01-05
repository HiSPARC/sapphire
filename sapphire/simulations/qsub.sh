#!/bin/bash -l
if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi
workon core
export JOB_HASH
python -m sapphire.simulations.qsub
