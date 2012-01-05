#!/bin/sh

if [ -d src ]; then
    #export PYTHONPATH=`pwd`/src:$PYTHONPATH
    export PYTHONPATH=`pwd`/src:/usr/local/lib/python
    #export PYTHONPATH=`pwd`/framework:$PYTHONPATH
fi
