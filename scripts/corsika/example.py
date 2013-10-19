#!/usr/bin/env python

import sys
import sapphire.corsika

if len(sys.argv) != 2 or sys.argv[1] == '-h':
    print "Usage: %s <filename>"%sys.argv[0]
    print 'For an example: %s $AUGEROFFLINEROOT/share/auger-offline/doc/SampleShowers/Corsika-1e19-6.part'%sys.argv[0]

    exit(0)

cors_file = corsika.CorsikaFile(sys.argv[1])

cors_file.check()

for event in cors_file.get_events():
    print event.get_header()
    print event.get_trailer()
    count = 0
    for particle in event.get_particles():
        count += 1

    print count, " particles"
