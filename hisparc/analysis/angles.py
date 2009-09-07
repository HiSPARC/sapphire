import zlib
from numpy import (mean, std, nan, sqrt, arctan2, cos, sin, arcsin, pi,
                   isnan, array, arccos, arcsin)
import pylab

ADC_THRESHOLD = 50
ADC_TIME_PER_SAMPLE = 2.5e-9

detector_positions = [(78.0, 15.05),
                      (78.0, 20.82),
                      (83.0, 23.71),
                      (73.0, 23.71)]

def reconstruct_angles(coincidence, dataset):
    """Reconstruct shower angles of a coincidence event"""

    times = reconstruct_time_differences(coincidence, dataset)
    return reconstruct_all_angles(times)

def reconstruct_all_angles(times):
    angles = []
    for detectors in xuniqueCombinations(range(4), 3):
        angles.append(reconstruct_angles_from_set(detectors, times))

    return angles

def reconstruct_time_differences(coincidence, dataset):
    """Reconstruct time differences of a coincidence event"""

    traces = dataset.root.hisparc.traces
    trace_idx = coincidence['hisparc_traces']

    times = []
    for i in range(4):
        trace = [int(x) for x in 
                    zlib.decompress(traces[trace_idx[i]]).split(',')[:-1]]
        times.append(reconstruct_time_from_trace(trace))

    return times

def reconstruct_time_from_trace(trace):
    """Reconstruct time of measurement from a trace"""

    t = trace[:100]
    baseline = mean(t)
    stdev = std(t)

    trace = trace - baseline
    threshold = ADC_THRESHOLD

    value = nan
    for i in range(len(trace)):
        if trace[i] >= threshold:
            value = i
            break

    return value * ADC_TIME_PER_SAMPLE

def xuniqueCombinations(items, n):
    """Return unique combinations from a list

    Reference: http://code.activestate.com/recipes/190465/

    """
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc

def reconstruct_angles_from_set(detectors, times):
    """Reconstruct angles from a set of three detectors"""

    c = 3.00e+8

    r1 = calc_distance(detectors[0], detectors[1])
    r2 = calc_distance(detectors[0], detectors[2])

    dt1 = times[detectors[0]] - times[detectors[1]]
    dt2 = times[detectors[0]] - times[detectors[2]]

    phi1 = calc_phi(detectors[0], detectors[1])
    phi2 = calc_phi(detectors[0], detectors[2])

    phi = arctan2((dt2 * r1 * cos(phi1) - dt1 * r2 * cos(phi2)),
                  (dt2 * r1 * sin(phi1) - dt1 * r2 * sin(phi2)) * -1)

    theta = arcsin(c * dt1 / (r1 * cos(phi - phi1)))

    if theta < 0:
        theta *= -1
        phi += pi
    
    phi %= 2 * pi

    #print detectors
    #print 'det 1:', dt1, phi1, r1
    #print 'det 2:', dt2, phi2, r2
    #print 'result:', theta, phi
    #t1 = r1 * cos(phi - phi1) * sin(theta) / 3e8
    #t2 = r2 * cos(phi - phi2) * sin(theta) / 3e8
    #print 'verification:', t1, t2

    #print dt2 * r1 * cos(phi1), dt1 * r2 * cos(phi2)
    #print dt2 * r1 * sin(phi1), dt1 * r2 * sin(phi2)
    
    return detectors, theta, phi

def calc_distance(det0, det1):
    """Calculate the distance between two detectors"""

    x0, y0 = detector_positions[det0]
    x1, y1 = detector_positions[det1]

    return sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)

def calc_phi(det0, det1):
    """Calculate the angle of the line of sight between two detectors"""

    x0, y0 = detector_positions[det0]
    x1, y1 = detector_positions[det1]

    return arctan2((y1 - y0), (x1 - x0))

def average_angles(angles):
    """Average a set of angles and estimate the error

    This function averages a set of angles using directional statistics and
    estimates the error using a custom function.

    Nan's are handled carefully. The returned estimate is in the interval
    [0, 1]. The closer it is to 1, the more accurate the average angle is.
    When it is zero, only Nan's have been averaged or the angles are
    separated by 180 degrees.

    """
    sines, cosines = [], []

    for angle in angles:
        if isnan(angle):
            sines.append(0)
            cosines.append(0)
        else:
            sines.append(sin(angle))
            cosines.append(cos(angle))

    mean_angle = arctan2(mean(sines), mean(cosines))
    radius = sqrt(mean(sines)**2 + mean(cosines)**2)

    return mean_angle, radius

def opening_angle(theta1, phi1, theta2, phi2):
    """Calculate the angle between two unit vectors"""

    return arccos(sin(theta1) * sin(theta2) * cos(phi1 - phi2) +
                  cos(theta1) * cos(theta2))

def angle_graphs(data, coincidences, ph_threshold, accuracy):
    """Produce some nice graphs"""

    cs = [x.fetch_all_fields() for x in coincidences if min(x['hisparc_pulseheights']) >= ph_threshold]
    angles = [reconstruct_angles(x, data) for x in cs]
    theta, phi = zip(*[[average_angles(y) for y in zip(*[(x[1], x[2]) for x in t])] for t in angles])

    ht, kt = [], []
    for i in range(len(theta)):
        if theta[i][1] > accuracy:
            ht.append(theta[i][0])
            kt.append(cs[i]['kascade_zenith'])
    
    hp, kp = [], []
    for i in range(len(phi)):
        if phi[i][1] > accuracy:
            hp.append(phi[i][0])
            kp.append(cs[i]['kascade_azimuth'])
    kp = array(kp)
    kp -= pi

    oa = []
    for i in range(len(angles)):
        if theta[i][1] > accuracy and phi[i][1] > accuracy:
            oa.append(opening_angle(theta[i][0], phi[i][0], cs[i]['kascade_zenith'], cs[i]['kascade_azimuth']))

    pylab.figure()
    pylab.hist(ht, bins=50, histtype='step', range=[0, pi/2], label='HiSPARC')
    pylab.hist(kt, bins=50, histtype='step', range=[0, pi/2], label='KASCADE')
    pylab.legend()
    pylab.xlabel('Zenith angle (radians)')
    pylab.ylabel('Count')
    pylab.title('Histogram of averaged zenith angles (accur > %.2f)' % accuracy)

    pylab.figure()
    pylab.hist(hp, bins=50, histtype='step', range=[-pi, pi], label='HiSPARC')
    pylab.hist(kp, bins=50, histtype='step', range=[-pi, pi], label='KASCADE')
    pylab.legend()
    pylab.xlabel('Azimuthal angle (radians)')
    pylab.ylabel('Count')
    pylab.title('Histogram of averaged azimuthal angles (accur > %.2f)' % accuracy)

    pylab.figure()
    pylab.hist(oa, bins=50, histtype='step')
    pylab.legend()
    pylab.xlabel('Angle (radians)')
    pylab.ylabel('Count')
    pylab.title('Histogram of opening angles (accur > %.2f)' % accuracy)

def test_angles():
    for theta in arange(0, 5) * (pi / 8):
        # phi = 0
        dt1 = test_calc_dt(5, theta)
        dt2 = test_calc_dt(10, theta)
        print "Theta, phi:", theta, 0.
        test_with_ran([dt1, dt1, 0., dt2])

        # phi = pi
        dt1 = test_calc_dt(5, theta)
        dt2 = test_calc_dt(10, theta)
        print "Theta, phi:", theta, pi
        test_with_ran([dt1, dt1, dt2, 0.])

        # phi = pi / 2
        dt1 = test_calc_dt(2.89, theta)
        dt2 = test_calc_dt(8.66, theta)
        print "Theta, phi:", theta, pi / 2
        test_with_ran([dt2, dt1, 0., 0.])

        # phi = 3 * pi / 2
        dt1 = test_calc_dt(5.77, theta)
        dt2 = test_calc_dt(8.66, theta)
        print "Theta, phi:", theta, 3 * pi / 2
        test_with_ran([0., dt1, dt2, dt2])

        print

def test_calc_dt(r, theta):
    return r * sin(theta) / 3e8

def test_with_ran(times):
    #times = [x + random() * 1e-11 for x in times]
    angles = reconstruct_all_angles(times)
    for angle in angles:
        print angle

def test_for_nan(angles):
    for theta, phi in angles:
        if isnan(theta):
            return True
    return False
