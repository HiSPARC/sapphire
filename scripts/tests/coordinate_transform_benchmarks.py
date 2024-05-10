"""
This is a benchmarking module developed for testing the coordinate
transformations by Ethan van Woerkom. The most important functions are
oldvsnew_diagram and pyephem_comp, respectively testing the old against the new
transformations and the new ones against pyephem.

transformspeeds tests the speed of the new transformations
"""

import datetime
import random
import time

import matplotlib.pyplot as plt
import numpy as np

from sapphire.transformations import celestial, clock
from sapphire.utils import angle_between


def transformspeeds():
    print('Running speeds for 100.000 transformations of the astropy functions:')
    a = np.array([(0, 0)] * 100000)

    t0 = time.clock()

    celestial.equatorial_to_horizontal_astropy(0, 0, 1_000_000_000, a)
    t1 = time.clock() - t0

    celestial.equatorial_to_zenithazimuth_astropy(0, 0, 1_000_000_000, a)
    t2 = time.clock() - t0

    celestial.zenithazimuth_to_equatorial_astropy(0, 0, 1_000_000_000, a)
    t3 = time.clock() - t0

    celestial.zenithazimuth_to_equatorial_astropy(0, 0, 1_000_000_000, a)
    t4 = time.clock() - t0

    print('EQ->HO, EQ-> ZA, HO->EQ, ZA->EQ runtimes:')

    print(t1, t2, t3, t4)


def angle_between_horizontal(azimuth1, altitude1, azimuth2, altitude2):
    """Calculate the angle between two horizontal coordinates

    Using the haversine formula,
    from: https://www.movable-type.co.uk/scripts/latlong.html

    :param azimuth#: Azimuth parts of the coordinates in radians.
    :param altitude#: Altitude parts of the coordinates in radians.
    :return: Angle between the two coordinates in radians.

    """
    zenith1, azimuth1 = celestial.horizontal_to_zenithazimuth(altitude1, azimuth1)
    zenith2, azimuth2 = celestial.horizontal_to_zenithazimuth(altitude2, azimuth2)
    dlat = zenith1 - zenith2
    dlon = azimuth2 - azimuth1
    a = np.sin(dlat / 2) ** 2 + np.sin(zenith1) * np.sin(zenith2) * np.sin(dlon / 2) ** 2
    angle = 2 * np.arcsin(np.sqrt(a))

    return angle


def oldvsnew_diagram():
    """
    Visual accuracy comparisons of old and new transformations.
    Compares the correlations between the transformations:
    equatorial_to_horizontal and equatorial_to_zenith_azimuth_astropy
    horizontal_to_equatorial and horizontal_to_zenith_azimuth_astropy
    Makes a histogram of the error differences
    between these two functions as well.
    The errors seem to be in the order of 1000 arcsec
    :return: None

    Ethan van Woerkom is responsible for the benchmarking functions;
    refer to him for when something is unclear
    """
    # make random frames, in correct angle range and from utc time 2000-2020
    # boxes for the four different transformation results
    etoha = []
    etoh = []
    htoe = []
    htoea = []
    straight = lambda x: x  # straight trendline function

    # Create the data sets for eq to az
    frames = [
        (
            random.uniform(-90, 90),
            random.uniform(-180, 180),
            random.randint(946684800, 1577836800),
            random.uniform(0, 2 * np.pi),
            random.uniform(-0.5 * np.pi, 0.5 * np.pi),
        )
        for _ in range(100)
    ]
    for i in frames:
        etoha.append(celestial.equatorial_to_zenithazimuth_astropy(i[0], i[1], i[2], [(i[3], i[4])])[0])
        etoh.append(celestial.equatorial_to_zenithazimuth(i[0], i[1], clock.utc_to_gps(i[2]), i[3], i[4]))
    # Data sets for hor to eq
    for i in frames:
        htoe.append(
            celestial.horizontal_to_equatorial(
                i[0],
                clock.utc_to_lst(datetime.datetime.utcfromtimestamp(i[2]), i[1]),
                i[4],
                i[3],
            ),
        )
        htoea.extend(celestial.horizontal_to_equatorial_astropy(i[0], i[1], i[2], [(i[3], i[4])]))

    # Make figs eq -> zenaz
    plt.figure(1)
    plt.suptitle('Zen/Az correlation in rads (equatorial_to_zenithazimuth)')

    zenrange = [0, np.pi]
    plt.subplot(211)
    plt.title('Zenith')
    plt.axis(zenrange * 2)
    plt.xlabel('New (Astropy)')
    plt.ylabel('Old')

    # Make figure and add 1:1 trendline

    plt.plot([co[0] for co in etoha], [co[0] for co in etoh], 'r.', zenrange, straight(zenrange), '-')

    plt.subplot(212)
    plt.title('Azimuth')
    azrange = [-np.pi, np.pi]
    plt.axis(azrange * 2)
    plt.xlabel('New (Astropy)')
    plt.ylabel('Old')
    # Make figure and add 1:1 trendline
    plt.plot([co[1] for co in etoha], [co[1] for co in etoh], 'b.', azrange, straight(azrange), '-')
    plt.tight_layout()  # Prevent titles merging
    plt.subplots_adjust(top=0.85)

    # Make histogram of differences
    plt.figure(2)
    # Take diff. and convert to arcsec
    nieuw = np.array(etoh) - np.array(etoha)
    nieuw *= 360 * 3600 / (2 * np.pi)

    plt.hist([i[0] for i in nieuw], bins=20)
    plt.title('Zenith Old-New Error (equatorial_to_zenithazimuth)')
    plt.xlabel('Error (arcsec)')
    plt.ylabel('Counts')

    plt.figure(3)
    plt.hist([i[1] for i in nieuw], bins=20)
    plt.title('Azimuth Old-New Error (equatorial_to_zenithazimuth)')
    plt.xlabel('Error (arcsec)')
    plt.ylabel('Counts')

    # Make histogram of differences using the absolute distance in arcsec
    # this graph has no wrapping issues
    plt.figure(7)
    nieuw = np.array([angle_between(etoh[i][0], etoh[i][1], etoha[i][0], etoha[i][1]) for i in range(len(etoh))])
    nieuw *= 360 * 3600 / (2 * np.pi)
    plt.hist(nieuw, bins=20)
    plt.title('ZEN+AZ Old-New Error (equatorial_to_zenithazimuth)')
    plt.xlabel('Error (arcsec)')
    plt.ylabel('Counts')

    # Make figs hor - > eq

    plt.figure(4)
    plt.suptitle('RA/DEC  correlation in rads (horizontal_to_equatorial)')
    altrange = [-0.5 * np.pi, 0.5 * np.pi]
    plt.subplot(211)
    plt.title('Declination')
    plt.axis(altrange * 2)
    plt.xlabel('New (Astropy)')
    plt.ylabel('Old')
    # Make figure and add 1:1 trendline
    plt.plot([co[1] for co in htoea], [co[1] for co in htoe], 'r.', altrange, straight(altrange), '-')

    plt.subplot(212)
    plt.title('Right Ascension')
    azrange = [0, 2 * np.pi]
    plt.axis(azrange * 2)
    plt.xlabel('New (Astropy)')
    plt.ylabel('Old')
    # Make figure and add 1:1 trendline
    plt.plot([co[0] for co in htoea], [co[0] for co in htoe], 'b.', azrange, straight(azrange), '-')
    plt.tight_layout()  # Prevent titles merging
    plt.subplots_adjust(top=0.85)

    # Make histogram of differences
    plt.figure(5)
    # Take diff. and convert to arcsec
    nieuw = np.array(htoe) - np.array(htoea)
    nieuw *= 360 * 3600 / (2 * np.pi)
    plt.hist([i[1] for i in nieuw], bins=20)
    plt.title('Declination Old-New Error (horizontal_to_equatorial)')
    plt.xlabel('Error (arcsec)')
    plt.ylabel('Counts')

    plt.figure(6)
    # Take diff. and convert to arcsec
    nieuw = np.array(htoe) - np.array(htoea)
    nieuw *= 360 * 3600 / (2 * np.pi)
    plt.hist([i[0] for i in nieuw], bins=20)
    plt.title('Right Ascension Old-New Error (horizontal_to_equatorial)')
    plt.xlabel('Error (arcsec)')
    plt.ylabel('Counts')

    # Make histogram of differences using the absolute distance in arcsec
    # this graph has no wrapping issues
    plt.figure(8)
    nieuw = np.array(
        [angle_between_horizontal(htoe[i][0], htoe[i][1], htoea[i][0], htoea[i][1]) for i in range(len(htoe))],
    )
    # Take diff. and convert to arcsec
    nieuw /= 2 / np.pi * 360 * 3600
    plt.hist(nieuw, bins=20)
    plt.title('RA+DEC Old-New Error (horizontal_to_equatorial)')
    plt.xlabel('Error (arcsec)')
    plt.ylabel('Counts')

    plt.show()


try:
    # This try-except block contains a pyephem accuracy benchmarking function.
    # It uses this structure to accommodate people without pyephem.
    import ephem

    def pyephem_comp():
        """
        This function compares the values from transformations done by our
        new astropy functions with the pyephem numbers. It generates
        correlation graphs between the new and old functions
        and histograms of the frequency of errors. Most errors do not much
        exceed 10 arcsec. There is complete correlation
        i.e. visually all points are on the same 1:1 line.
        These comparisons are done on the basis of 100 randomly generated
        points
        :return: None

        Ethan van Woerkom is responsible for the benchmarking functions;
        refer to him for when something is unclear
        """
        # Set up randoms equatorial J2000 bodies
        # that we will convert the RAs/Decs of.
        eq = []  # random frames to use
        for i in range(100):
            eq.append(
                (
                    random.uniform(-90, 90),
                    random.uniform(-180, 180),
                    random.randint(946684800, 1577836800),
                    random.uniform(0, 2 * np.pi),
                    random.uniform(-0.5 * np.pi, 0.5 * np.pi),
                ),
            )
        efemeq = []  # store pyephem transformations to equatorial
        altaz = []  # store pyephem transformations to altaz (horizontal)
        htoea = []  # store astropy transformations to equatorial
        etoha = []  # store astropy transformations to horizontal (altaz)
        for latitude, longitude, utc, ra, dec in eq:
            # Calculate altaz
            # Set observer for each case
            obs = ephem.Observer()
            obs.lat = str(latitude)
            obs.lon = str(longitude)
            obs.date = datetime.datetime.utcfromtimestamp(utc)
            obs.pressure = 0  # Crucial to prevent refraction correction!

            # Set body for each case
            coord = ephem.FixedBody()
            coord._ra = ra
            coord._dec = dec

            # Do calculation and add to list, it is necessary
            # to force the coords into radians
            coord.compute(obs)
            altaz.append((float(coord.az), float(coord.alt)))

            # Also calculate efemeq using eq
            result = obs.radec_of(ra, dec)  # This is of course not ra,dec but
            # actually az, alt.

            efemeq.append((float(result[0]), float(result[1])))

        # Produce horizontal_to_equatorial_astropy results
        for latitude, longitude, utc, az, alt in eq:
            result = celestial.horizontal_to_equatorial_astropy(latitude, longitude, utc, [(az, alt)])
            htoea.extend(result)

        # Produce equatorial_to_horizontal_astropy results
        for latitude, longitude, utc, ra, dec in eq:
            result = celestial.equatorial_to_horizontal_astropy(latitude, longitude, utc, [(ra, dec)])
            etoha.extend(result)

        altdecrange = [-0.5 * np.pi, 0.5 * np.pi]
        azrarange = [0, 2 * np.pi]

        plt.figure(1)
        plt.suptitle('RA/Dec correlation Altaz->(Astropy/Pyephem)->RA,DEC')

        # Create RA correlation subplot
        plt.subplot(211)
        plt.title('RA')
        plt.axis(azrarange * 2)
        plt.xlabel('Pyephem RA (rad)')
        plt.ylabel('Astropy RA (rad)')

        # Plot RA correlation and trendline
        plt.plot([co[0] for co in efemeq], [co[0] for co in htoea], 'b.', azrarange, azrarange, '-')

        # DEC correlation subplot
        plt.subplot(212)
        plt.title('DEC')
        plt.axis(altdecrange * 2)
        plt.xlabel('Pyephem DEC (rad)')
        plt.ylabel('Astropy DEC (rad)')

        # Plot DEC correlation and trendline
        plt.plot([co[1] for co in efemeq], [co[1] for co in htoea], 'r.', altdecrange, altdecrange, '-')

        # Formatting
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)

        # Make RA error histogram
        plt.figure(2)
        plt.title('RA Error Altaz->(Astropy/Pyephem)->RA,DEC')

        nieuw = np.array(htoea) - np.array(efemeq)
        # Get differences in arcsec
        nieuw *= 360 * 3600 / (2 * np.pi)

        plt.hist([co[0] for co in nieuw], bins=20)

        plt.ylabel('Counts')
        plt.xlabel('Error (arcsec)')

        # Make DEC error histogram
        plt.figure(3)
        plt.title('DEC Error Altaz->(Astropy/Pyephem)->RA,DEC')
        plt.hist([co[1] for co in nieuw], bins=20)

        plt.ylabel('Counts')
        plt.xlabel('Error (arcsec)')

        # Altaz comparison plot
        plt.figure(4)
        plt.suptitle('Alt/Az correlation RA,DEC->(pyephem/astropy)->Altaz')

        # Altitude
        plt.subplot(211)
        plt.title('Altitude')
        plt.axis(altdecrange * 2)
        plt.xlabel('Pyephem Altitude (rad)')
        plt.ylabel('Astropy Altitude (rad')

        # Plot with trendline
        plt.plot([co[1] for co in altaz], [co[1] for co in etoha], 'b.', altdecrange, altdecrange, '-')

        # Azimuth
        plt.subplot(212)
        plt.title('Azimuth')
        plt.axis(azrarange * 2)
        plt.xlabel('Pyephem Azimuth (rad)')
        plt.ylabel('Astropy Azimuth (rad)')

        plt.plot([co[0] for co in altaz], [co[0] for co in etoha], 'r.', azrarange, azrarange, '-')

        # Formatting
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)

        # Alt error histogram
        plt.figure(5)
        plt.title('Altitude Error RA,DEC->(pyephem/astropy)->Altaz')
        nieuw = np.array(etoha) - np.array(altaz)
        nieuw *= 360 * 3600 / (2 * np.pi)
        plt.hist([co[1] for co in nieuw], bins=20)

        plt.ylabel('Counts')
        plt.xlabel('Error (arcsec)')

        # Az error histogram
        plt.figure(6)
        plt.title('Azimuth Error RA,DEC->(pyephem/astropy)->Altaz')
        plt.hist([co[0] for co in nieuw], bins=20)

        plt.ylabel('Counts')
        plt.xlabel('Error (arcsec)')

        # Make histograms of differences using the absolute distance in arcsec
        # these graphs have no wrapping issues

        plt.figure(7)
        nieuw = np.array(
            [angle_between_horizontal(altaz[i][0], altaz[i][1], etoha[i][0], etoha[i][1]) for i in range(len(etoha))],
        )
        nieuw *= 360 * 3600 / (2 * np.pi)
        plt.hist(nieuw, bins=20)
        plt.title('Alt+Azi Error RA,DEC->(pyephem/astropy)->Altaz')
        plt.xlabel('Error (arcsec)')
        plt.ylabel('Counts')

        plt.figure(8)
        nieuw = np.array(
            [angle_between_horizontal(efemeq[i][0], efemeq[i][1], htoea[i][0], htoea[i][1]) for i in range(len(htoea))],
        )
        # Take difference and convert to arcsec
        nieuw *= 360 * 3600 / (2 * np.pi)

        plt.hist(nieuw, bins=20)
        plt.title('RA+DEC Error Altaz->(pyephem/astropy)->RA,DEC')
        plt.xlabel('Error (arcsec)')
        plt.ylabel('Counts')

        # Done; output
        plt.show()

except ImportError:
    # Pyephem is not required so there is a case for when it is not present
    def pyephem_comp():
        print('Pyephem not present; no comparisons will be done')
