import tables
import csv
from itertools import izip

try:
    data
except NameError:
    data = tables.openFile('kascade.h5', 'a')

events = data.root.efficiency.events

try:
    f.close()
except NameError:
    pass
finally:
    f = open('swaarde+.csv', 'r')

reader = csv.reader(f, delimiter='\t')

for table_row, csv_row in izip(events, reader):
    Ne, Ze, x, y, d0, d1, d2, d3, s = [float(u) for u in csv_row]

    assert Ne == table_row['k_Num_e']
    assert Ze == table_row['k_zenith']
    assert ((x, y) == table_row['k_core_pos']).all()
    #assert ((d0, d1, d2, d3) == table_row['k_dens_e']).all()

    table_row['s'] = s
    table_row.update()

events.flush()
