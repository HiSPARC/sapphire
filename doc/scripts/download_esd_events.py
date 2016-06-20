import datetime
import tables
from sapphire import esd

DATAFILE = 'data.h5'
STATIONS = [501, 503, 506]
START = datetime.datetime(2016, 1, 1)
END = datetime.datetime(2016, 1, 2)


if __name__ == '__main__':
    if 'data' not in globals():
        data = tables.open_file(DATAFILE, 'a')

    for station in STATIONS:
        group = '/s%d' % station
        if group not in data:
            esd.download_data(data, group, station, START, END)
