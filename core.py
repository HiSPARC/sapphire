import subprocess
from math import sqrt

p = subprocess.Popen(['gunzip', '-c', 'HiSparc.dat.gz'],
                     stdout=subprocess.PIPE)

x, y, r = [], [], []

i = 0
while True:
    line = p.stdout.readline()
    if not line:
        break

    data = line.split(' ')
    data = [float(v) for v in data]

    Irun, Ieve, Gt, Mmn, EnergyArray, Xc, Yc, Ze, Az, Size, Nmu, He0, \
    Hmu0, He1, Hmu1, He2, Hmu2, He3, Hmu3, P200, T200 = data

    #if Gt >= 1214918985:
    x.append(Xc)
    y.append(Yc)
    r.append(sqrt((Xc-80)**2+(Yc-20)**2))
    i += 1

    if i > 1000:
        break
