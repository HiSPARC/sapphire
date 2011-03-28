import tables
import progressbar as pb
import csv
data = tables.openFile('kascade.h5', 'r')
events = data.root.efficiency.events
f = open('swaarde.csv', 'w')
writer = csv.writer(f, delimiter='\t')
progress = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(), pb.ETA()])
for e in progress(events):
    Ne = e['k_Num_e']
    Ze = e['k_zenith']
    x, y = e['k_core_pos']
    d0, d1, d2, d3 = e['k_dens_e']
    writer.writerow((Ne, Ze, x, y, d0, d1, d2, d3))
