import tables

from hisparc.containers import KascadeEvent


DATAFILE = 'kascade.h5'
GROUP = '/efficiency'

DETECTORS = [(65., 15.05, 'UD'), (65., 20.82, 'UD'), (70., 23.71, 'LR'),
             (60., 23.71, 'LR')]
LOW_THRES = round(30 / .57)
HIGH_THRES = round(70 / .57)
TRIGGER = lambda l, h: l >= 3 or h >= 2

Event = {'core_dist': tables.FloatCol(),
         'core_alpha': tables.FloatCol()}
Trigger = {'n_low': tables.UInt8Col(),
           'n_high': tables.UInt8Col(),
           'self_triggered': tables.BoolCol()}

def main(data):
    if not GROUP in data:
        data.createGroup('/', GROUP[1:])

    copy_and_add_to_table(data, '/kascade/events', GROUP, 'kascade')
    copy_and_add_to_table(data, '/coincidences/events', GROUP, 'events')

def copy_and_add_to_table(data, src_node, dst_group, dst_node):
    src_node = data.getNode(src_node)
    dst_group = data.getNode(dst_group)

    if dst_node in dst_group:
        data.removeNode(dst_group, dst_node)

    table_desc = src_node.coldescrs.copy()
    table_desc.update(Event)
    if 'pulseheights' in src_node.colnames:
        table_desc.update(Trigger)
        analyze_trigger = True
    else:
        analyze_trigger = False
    table = data.createTable(dst_group, dst_node, table_desc)

    dst_row = table.row
    if 'core_pos' in src_node.colnames:
        core_pos = 'core_pos'
    else:
        core_pos = 'k_core_pos'
    for row in src_node:
        for col in src_node.colnames:
            dst_row[col] = row[col]
        x, y = row[core_pos]
        x -= DETECTORS[1][0]
        y -= DETECTORS[1][1]
        dst_row['core_dist'] = sqrt(x ** 2 + y ** 2)
        dst_row['core_alpha'] = arctan2(y, x)
        if analyze_trigger:
            n_low = (row['pulseheights'] >= LOW_THRES).sum()
            n_high = (row['pulseheights'] >= HIGH_THRES).sum()
            dst_row['n_low'] = n_low
            dst_row['n_high'] = n_high
            dst_row['self_triggered'] = TRIGGER(n_low, n_high)
        dst_row.append()
    table.flush()
    data.flush()


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('kascade.h5', 'a')

    main(data)
