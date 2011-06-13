#!/bin/sh
if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi
Aires < sime15.inp
