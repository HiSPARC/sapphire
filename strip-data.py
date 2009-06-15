import gzip

if __name__ == '__main__':
    g = gzip.open('HiSparc.dat.gz')
    f = open('kascade.dat', 'w')

    while True:
        # read a line from the subprocess stdout buffer
        line = g.readline()
        if not line:
            # no more lines left, EOF
            break

        # break up the line into an array of floats
        data = line.split(' ')
        data = [float(x) for x in data]

        # read all columns into KASCADE-named variables
        Irun, Ieve, Gt, Mmn, EnergyArray, Xc, Yc, Ze, Az, Size, Nmu, He0, \
        Hmu0, He1, Hmu1, He2, Hmu2, He3, Hmu3, P200, T200 = data

        if Gt < 1214784000:
            # Before June 30th, 2008
            continue
        elif Gt > 1215043200:
            # After July 2th, 2008
            break
        else:
            f.write(line)
