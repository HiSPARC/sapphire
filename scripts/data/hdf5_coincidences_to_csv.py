#!/usr/bin/env python

"""Convert a HDF5 file with coincidences to csv

Use this script to convert a HDF5 data file with events and coincidences
to csv. The output will contain a coincidence row followed by N (numbe
rof events in the coincidence) event rows.

The output will contain two header rows, first the header for the
coincidence rows, followed by the header for the event rows.

Currently this only works for data files where coincidences have been
stored by :class:`sapphire.analysis.coincidences.CoincidencesESD`

"""
import csv

import tables


def data_to_csv(destination, source, coincidences_group):
    """Copy data from the source HDF5 file to the destination csv

    :param destination: Writable CSV file
    :param source: Readable PyTables file
    :param coincidences_group: Path in the HDF5 file to coincidences group

    """
    csvwriter = csv.writer(destination, delimiter='\t')
    c_group = source.getNode(coincidences_group)
    coincidences = c_group.coincidences
    c_index = c_group.c_index
    s_index = c_group.s_index

    try:
        cluster = c_group._v_attrs.cluster
        s_numbers = [station.number for station in cluster.stations]
    except:
        s_numbers = [s_group.split('_')[-1] for s_group in s_index[:]]

    csvwriter.writerow(coincidences.colnames)
    csvwriter.writerow(['station_number'] +
                       source.getNode(s_index[0]).events.colnames)

    for coincidence in coincidences[:]:
        csvwriter.writerow(coincidence)
        c_idx = c_index[coincidence['id']]
        for s_idx, e_idx in c_idx:
            event = source.getNode(s_index[s_idx], 'events')[e_idx]
            csvwriter.writerow([s_numbers[s_idx]] + list(event))


if __name__ == "__main__":
    # Path to the csv file to be created.
    destination = 'data.csv'

    # Path to the HDF5 file to read.
    source = 'data.h5'

    # Path in the HDF5 to the group where the coincidences, c_index and
    # s_index tables are located.
    coincidences_group = '/network/coincidences'

    with open(destination, 'wb') as csvfile:
        with tables.openFile(source, 'r') as hdf5data:
            data_to_csv(csvfile, hdf5data, coincidences_group)
