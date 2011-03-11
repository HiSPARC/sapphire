import tables
from numpy import sqrt, arctan2, seterr, deg2rad, pi, isnan, cos
import progressbar as pb

from hisparc.containers import KascadeEvent
from analysis_kascade import reconstruct_angle, ADC_MIP


DATAFILE = 'kascade.h5'
GROUP = '/efficiency'

DETECTORS = [(65., 15.05, 'UD'), (65., 20.82, 'UD'), (70., 23.71, 'LR'),
             (60., 23.71, 'LR')]
LOW_THRES = round(30 / .57)
HIGH_THRES = round(70 / .57)
TRIGGER = lambda l, h: l >= 3 or h >= 2

Event = {'core_dist': tables.FloatCol(),
         'core_alpha': tables.FloatCol(),
         'k_cosdens_e': tables.FloatCol(shape=4)}
Trigger = {'n_low': tables.UInt8Col(),
           'n_high': tables.UInt8Col(),
           'self_triggered': tables.BoolCol()}
Reconstruction = {'n': tables.FloatCol(shape=4),
                  't': tables.FloatCol(shape=4),
                  'h_theta': tables.FloatCol(),
                  'h_phi': tables.FloatCol(),
                  'k_theta': tables.FloatCol(),
                  'k_phi': tables.FloatCol(),
                  'reconstructed': tables.BoolCol()}


def main(data):
    seterr(invalid='ignore')
    if not GROUP in data:
        data.createGroup('/', GROUP[1:])

    build_efficiency_dataset(data, '/coincidences/events', GROUP, 'events')

def build_efficiency_dataset(data, src_node, dst_group, dst_node):
    src_node = data.getNode(src_node)
    dst_group = data.getNode(dst_group)

    timing_data = data.root.analysis.timing_data.read()

    if dst_node in dst_group:
        data.removeNode(dst_group, dst_node)

    table_desc = src_node.coldescrs.copy()
    table_desc.update(Event)
    table_desc.update(Trigger)
    table_desc.update(Reconstruction)
    table = data.createTable(dst_group, dst_node, table_desc)

    dst_row = table.row
    progress = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(),
                                       pb.ETA()])
    for row, timing in progress(zip(src_node[:], timing_data[:])):
        for col in src_node.colnames:
            dst_row[col] = row[col]
        x, y = row['k_core_pos']
        x -= DETECTORS[1][0]
        y -= DETECTORS[1][1]
        dst_row['core_dist'] = sqrt(x ** 2 + y ** 2)
        dst_row['core_alpha'] = arctan2(y, x)
        n_low = (row['pulseheights'] >= LOW_THRES).sum()
        n_high = (row['pulseheights'] >= HIGH_THRES).sum()
        dst_row['n_low'] = n_low
        dst_row['n_high'] = n_high
        dst_row['self_triggered'] = TRIGGER(n_low, n_high)

        n = row['pulseheights'] / ADC_MIP

        event = dict(n1=n[0], n2=n[1], n3=n[2], n4=n[3], t1=timing[0],
                     t2=timing[1], t3=timing[2], t4=timing[3])
        theta, phi, theta1, theta2 = reconstruct_angle(event)

        dst_row['k_cosdens_e'] = row['k_dens_e'] * cos(row['k_zenith'])
        dst_row['k_theta'] = row['k_zenith']
        dst_row['k_phi'] = -(row['k_azimuth'] + deg2rad(75)) % \
                           (2 * pi) - pi

        dst_row['n'] = n
        dst_row['t'] = timing
        dst_row['h_theta'] = theta
        dst_row['h_phi'] = phi
        if not isnan(theta) and not isnan(phi):
            dst_row['reconstructed'] = True
        else:
            dst_row['reconstructed'] = False

        dst_row.append()
    table.flush()
    data.flush()


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('kascade.h5', 'a')

    main(data)
