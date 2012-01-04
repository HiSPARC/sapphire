#!/bin/sh

if [ -d src ]; then
    #export PYTHONPATH=`pwd`/src:$PYTHONPATH
    export PYTHONPATH=`pwd`/src
    export PYTHONPATH=`pwd`/framework:$PYTHONPATH
fi
