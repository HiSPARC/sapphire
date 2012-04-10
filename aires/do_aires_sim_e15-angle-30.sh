#!/bin/bash -l
if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi
Aires < sime15-angle-30.inp
