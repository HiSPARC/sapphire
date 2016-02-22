import datetime
import tables
from sapphire import esd
from sapphire.analysis import reconstructions

DATAFILE = 'data.h5'
STATIONS = [501, 503, 506]
START = datetime.datetime(2016, 1, 1)
END = datetime.datetime(2016, 1, 2)


if __name__ == '__main__':
    if 'data' not in globals():
        data = tables.open_file(DATAFILE, 'a')

    if '/coincidences' not in data:
        esd.download_coincidences(data, stations=STATIONS, start=START, end=END)

    if '/reconstructions' not in data:
        rec = reconstructions.ReconstructESDCoincidences(data)
        rec.reconstruct_and_store()
