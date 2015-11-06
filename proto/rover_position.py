#!/usr/bin/python
#! encoding: UTF8

from parse_rinex import parse_rinex, utc2gpst, gps2utc, LEAP, A0, A1
import numpy as np
from numpy import sin, cos, sqrt

__author__ = 'kirienko'

"""
1. Parse Nav file --> ephemerids --> sat coords (ECEF)
2. Parse Obs file --> pseudoranges to rover
3. Calculate rover's coords in ECEF
"""

# Pretty printing
def lla_string(R):
    return "φ = %.2f°, θ = %.2f°, h = %d m" % (R[0],R[1],int(R[2]))
def xyz_string(R):
    return "(%d, %d, %d) [km]" % tuple(map(int,R / 1000))

def ecef_to_lat_lon_alt1(R):
    """
    Fukushima implementation of the Bowring algorithm,
    see [3] -- equations (C7) - (C12)
    :param R: (X,Y,Z) -- coordinates in ECEF (numpy array)
    :return:  (φ,θ,h) -- lat [deg], lon [deg], alt in WGS84 (numpy array)
    """
    # WGS 84 constants
    a = 6378137.0       # Equatorial Radius [m]
    b = 6356752.314245  # Polar Radius [m]
    e_sq = 0.08181919092890624 # e = sqrt(1-b²/a²)
    # e1 = b/a            # e' = sqrt(1 - e²)
    # print e1, sqrt(1-e_sq)
    e1 = sqrt(1-e_sq)            # e' = sqrt(1 - e²)
    c = a * e_sq        # (6)
    if isinstance(R,list): R = np.array(R)
    p = sqrt(R[0]**2 + R[1]**2)             # (1) - sqrt(X² + Y²)
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
        h = sqrt(T1 - e_sq)/e1 * (p - a/sqrt(T1))
        print "p > z, p = %d, z = %d, a/sqrt(T1) =" % (int(p),int(R[2])), a/sqrt(T1)
    else:                                   # (C12)
        h = sqrt(T1 - e_sq) * (R[2]/T - b/sqrt(T1))
        print "p < z"
    out = np.array([np.degrees(phi), np.degrees(theta), h])
    print "DEBUG: %s --> %s" % (xyz_string(R),lla_string(out))
    # return np.array([theta, phi, h])
    # return
    return out

def ecef_to_lat_lon_alt(R):
    """
    Fukushima implementation of the Bowring algorithm (2006),
    see [?] --
    :param R: (X,Y,Z) -- coordinates in ECEF (numpy array)
    :return:  (φ,θ,h) -- lat [deg], lon [deg], alt in WGS84 (numpy array)
    """
    # WGS 84 constants
    a = 6378137.0       # Equatorial Radius [m]
    # b = 6356752.314245  # Polar Radius [m]
    E = 0.08181919092890624 # e = sqrt(1-b²/a²)
    e1 = sqrt(1 - E)            # e' = sqrt(1 - e²)
    b = a * e1
    if isinstance(R,list): R = np.array(R)
    p = sqrt(R[0]**2 + R[1]**2)             # (1) - sqrt(X² + Y²)
    az = abs(R[2])
    Z = e1*az/a
    P = p/a
    S, C = Z, e1 * P             # (C8) - zero approximation
    max_iter = 5
    for i in range(max_iter):
        A = sqrt(S**2 + C**2)
        Cn = P * A**3 - E * C**3
        Sn = Z * A**3 + E * S**3
        delta = abs(Sn/Cn - S/C)*C/S
        if abs(delta) < 1e-10 or i == max_iter-1:
            break
        S, C = Sn, Cn
    # print "i =",i
    theta = np.math.atan2(R[1],R[0])
    Cc = e1 * Cn
    phi = np.sign(R[2])*np.math.atan2(Sn,Cc)
    h = (p*Cc + az*Sn - b*sqrt(Sn**2 + Cn**2))/sqrt(Cc**2 + Sn**2)
    out = np.array([np.degrees(phi), np.degrees(theta), h])
    # print "\nDEBUG: %s --> %s" % (xyz_string(R),lla_string(out))
    # return np.array([theta, phi, h])
    return out


def lat_lon_alt_to_ecef_xyz(R):
    """
    see [1] p. 512
    :param R: (φ,θ,h) -- lat [deg], lon [deg], alt in WGS84 (numpy array)
    :return:  (X,Y,Z) -- coordinates in ECEF (numpy array)
    """
    # WGS 84 constants
    a2 = 6378137.0**2       # Equatorial Radius [m]
    b2 = 6356752.314245**2  # Polar Radius [m]
    radius = lambda phi: sqrt(a2*(cos(phi))**2 + b2*(sin(phi))**2)
    f,t,h = np.deg2rad(R[0]),np.deg2rad(R[1]),R[2]
    X = cos(t)*cos(f)*(h + a2/radius(f))
    Y = sin(t)*cos(f)*(h + a2/radius(f))
    Z = sin(f)*(h + b2/radius(f))
    return np.array([X,Y,Z])


def least_squares(obs, navs):
    """
    x = (A^TA)^{-1}A^T l
    Takes observation
    :return:
    """
    c = 299792428 # speed of light
    print "=== Least Squares ==="
    # from A-matrix
    # ↳ find all possible satellites N
    sats = []
    for i,r in enumerate(obs.PRN_number):
        if obs.obs_data['C1'][i] and obs.obs_data['P2'][i] and 'G' in r:
            sats += [r]
    # ↳ if N >= 4: form matrix
    if len(sats) > 3:
        # observed [iono-free] pseudoranges
        P = np.array([obs.ionofree_pseudorange(s) for s in sats])
        # print "P =",P
        # get XYZ-coords of satellites
        XYZs = np.array([navs[s][0].eph2pos(utc2gpst(o.date)) for s in sats])
        # print "XYZs =",XYZs

    xyz = np.zeros(4) # initial point
    for itr in range(20):
        # geometrical ranges
        rho = np.array([np.sqrt(sum([(x - xyz[i])**2 for i,x in enumerate(XYZs[j])])) for j in range(len(sats))])
        # print "rho =",rho

        # from A-matrix
        A = np.matrix([np.append((xyz[:3] - XYZs[i])/rho[i],1) for i in range(len(sats))])
        # print "A:\n",A
        AT = A.transpose()
        # form l-vector
        l = np.matrix([P[i] - rho[i] + xyz[3] for i in xrange(len(sats))]).transpose()
        # form x-vector
        x_hat = ((AT*A).I * AT * l).flatten().getA()[0]
        # x_hat[3] /= c
        print "(x,y,z,cδt) =",x_hat
        # iterate
        xyz += x_hat
        delta = np.sqrt(sum(map(lambda k: k**2,x_hat[:3])+[x_hat[3]/c]))
        if delta < 10.:
            break
    return xyz

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
    o = observations[240]
    print o.sat_types
    # sats = ['G05', 'G16', 'G18', 'G21']
    # for s in sats:
    #     print "C1 from %s at time %s: %f [m]"% (s,str(o.date.time()),o.pseudorange(s,'C1'))


    print lla_string(ecef_to_lat_lon_alt(least_squares(o,navigations)))

    # ecef_to_lat_lon_alt([3695041.,3695041.,3695041.])   # 45,45,0
    # ecef_to_lat_lon_alt([6400000*.707,6400000*.707,1])  # 0.45,0
    # ecef_to_lat_lon_alt([1,1,6400000])
    # ecef_to_lat_lon_alt([1,6400000,1])
    # ecef_to_lat_lon_alt([6400000,1,1])