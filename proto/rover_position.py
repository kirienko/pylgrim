#!/usr/bin/python
#! encoding: UTF8

from parse_rinex import parse_rinex, utc2gpst, gps2utc, LEAP, A0, A1
import numpy as np

__author__ = 'kirienko'

"""
1. Parse Nav file --> ephemerids --> sat coords (ECEF)
2. Parse Obs file --> pseudoranges to rover
3. Calculate rover's coords in ECEF
"""



def ecef_to_lot_lon_alt(R): # TODO -- make iterative
    """
    Fukushima implementation of the Bowring algorithm,
    see [3] -- equations (C7) - (C12)
    :param R: (X,Y,Z) -- coordinates in ECEF (numpy array)
    :return:  (θ,φ,h) -- lat, lon, alt in WGS84 (numpy array)
    """
    # WGS 84 constants
    a = 6378137.0       # Equatorial Radius [m]
    b = 6356752.314245  # Polar Radius [m]
    e_sq = 0.08181919092890624 # e = sqrt(1-b²/a²)
    e1 = b/a            # e' = sqrt(1 - e²)
    c = a*e_sq          # a² - b²

    p = np.sqrt(np.square(R[:2]).sum())     # (1) - sqrt(X² + Y²)
    T = R[2]/(e1*p)              # (C8) - zero approximation
    for i in range(5):
        C = np.power(1 + T**2,-0.5)         # (C9)
        S = C * T                           # (C9)
        T_new = (e1 * R[2] + c * S**3)/(p - c * C**3)   # (C7)
        delta = T_new - T
        T = T_new
        if abs(delta)/T < 1e-9:
            break
    theta = np.math.atan2(R[1],R[0])
    phi = np.math.atan2(T,e1)               # (C10)
    T1 = 1 + T**2
    if p > R[2]:                            # (C11)
        h = np.sqrt(T1 - e_sq)/e1 * (p - a/np.sqrt(T1))
    else:                                   # (C12)
        h = np.sqrt(T1 - e_sq) * (R[2]/T - b/np.sqrt(T1))
    # print "%.2f, %.2f, %d m" % (np.degrees(theta),np.degrees(phi),int(h))
    return np.array([theta, phi, h])
    # return np.array([np.degrees(theta), np.degrees(phi), h])

def dummy_rover_pos(rho,R):
    """
    Naïve approach explained in [1], 2.4 (p. 48):
    calculate rover's position using measurement of 4 pseudoranges to 4 satellites.
    :param: rho = (r1,r2,r3,r4) -- pseudoranges to sat-s 1,2,3,4
    :param: R = (R1,R2,R3,R4) -- coordinates  of sat-s 1,2,3,4 at the moment of measurement
            Ri = (Xi,Yi,Zi) -- coordinates of i-th satellite in ECEF
    :return: X = M⁻¹ × Y
    """
    a = 6378137.0       # Equatorial Radius [m]
    b = 6356752.314245  # Polar Radius [m]
    ear = (a+b)/2  # earth's radius #FIXME
    print "rho", rho
    print "R0",R
    sq = lambda x: np.square(x).sum()
    Y = np.matrix([rho[i]-sq(R[i])-ear for i in range(4)]).transpose()
    M = np.matrix([np.append(-2*R[i],1) for i in range(4)])
    return (M.I*Y).transpose()

if __name__ == "__main__":

    nav_file = '../test_data/test.n'
    obs_file = '../test_data/test.o'

    # Process Nav file:
    # ``navigations`` is a dict with
    #       keys:   GNSS identificators, i.e. 'G16', 'G04', ...
    #       values: Nav observation objects
    navigations = parse_rinex(nav_file)

    # Process Obs file
    observations = parse_rinex(obs_file)
    o = observations[250]
    print o.sat_types
    sats = ['G05', 'G16', 'G18', 'G21']
    for s in sats:
        print "C1 from %s at time %s: %f [m]"% (s,str(o.date.time()),o.pseudorange(s,'C1'))

    Rho = np.array([o.pseudorange(s,'C1') for s in sats])
    # for s in sats:
    #     xyz = navigations[s][0].eph2pos(utc2gpst(o.date))
    #     print xyz
    XYZs = np.array([navigations[s][0].eph2pos(utc2gpst(o.date)) for s in sats])
    coords = (dummy_rover_pos(Rho,XYZs)).tolist()[0]
    t,f,h = ecef_to_lot_lon_alt(coords[:3])
    print "rover here: %.2f, %.2f, %d " % (np.degrees(t), np.degrees(f),h)