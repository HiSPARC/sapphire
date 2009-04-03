import subprocess
from math import sqrt

p = subprocess.Popen(['gunzip', '-c', 'HiSparc.dat.gz'],
                     stdout=subprocess.PIPE)
f = open('kascade-time.dat', 'w')

i = 0
while True:
    line = p.stdout.readline()
    if not line:
        break

    data = line.split(' ')
    data = [float(x) for x in data]

    Irun, Ieve, Gt, Mmn, EnergyArray, Xc, Yc, Ze, Az, Size, Nmu, He0, \
    Hmu0, He1, Hmu1, He2, Hmu2, He3, Hmu3, P200, T200 = data

    if Gt >= 1214918985:
        f.write("%d %d %e %f %f\n" % (Gt, Mmn, EnergyArray, He0, Hmu0))

f.close()
