#!/bin/sh

if [ -d src ]; then
    export PYTHONPATH=`pwd`/src:$PYTHONPATH
fi
