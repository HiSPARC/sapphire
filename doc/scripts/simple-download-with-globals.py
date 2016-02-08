import datetime
import tables
from sapphire import esd

DATAFILE = 'data.h5'
START = datetime.datetime(2016, 1, 1)
END = datetime.datetime(2016, 1, 2)


if 'data' not in globals():
    data = tables.open_file(DATAFILE, 'w')
    esd.download_data(data, '/s501', 501, START, END)
