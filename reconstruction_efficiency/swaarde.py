import tables
import progressbar as pb
import csv
data = tables.openFile('kascade.h5', 'r')
events = data.root.efficiency.events
f = open('swaarde.csv', 'w')
writer = csv.writer(f, 'excel-tab')
progress = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(), pb.ETA()])
for e in progress(events):
    Ne = e['k_Num_e']
    Ze = e['k_zenith']
    x, y = e['k_core_pos']
    d0, d1, d2, d3 = e['k_dens_e']
    p0, p1, p2, p3 = e['pulseheights']
    trig = e['self_triggered']
    writer.writerow((Ne, Ze, x, y, d0, d1, d2, d3, p0, p1, p2, p3, trig))
