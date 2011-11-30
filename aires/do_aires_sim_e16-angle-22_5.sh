#!/bin/bash -l
if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi
Aires < sime16-angle-22_5.inp
