#!/usr/bin/env python

import sys
from sapphire import corsika


cors_file = corsika.CorsikaFile('../../sapphire/tests/corsika/DAT000000')

cors_file.check()

for event in cors_file.get_events():
    print event.get_header()
    print event.get_end()
    count = 0
    for particle in event.get_particles():
        count += 1

    print count, " particles"
