""" HiSPARC Framework

    This collection of scripts intends to be useful for all aspects of
    HiSPARC data analysis.

    At the moment, it is nothing like a *real* framework, but hopefully,
    someday...

    The following packages and modules are included:

    :mod:`~hisparc.publicdb`
        a module to fetch data from the datastore

    :mod:`~hisparc.analysis`
        all modules related to reconstruction and analysis

    :mod:`~hisparc.gpstime`
        a module to convert GPS to UTC times and vice versa

    :mod:`~hisparc.kascade`
        a module to get data from KASCADE and search for coincidences
        between HiSPARC and KASCADE

"""
import containers
import publicdb
import kascade
import analysis
