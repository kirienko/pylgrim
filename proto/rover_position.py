#!/usr/bin/python
#! encoding: UTF8

import datetime as dt
import numpy as np

from coord.ecef import ecef_to_lat_lon_alt, sat_elev
from parse_rinex import parse_rinex, utc2gpst
from visualization.ellipsoid import satellites

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


def least_squares(obs, navs):
    """
    x = (A^TA)^{-1}A^T l
    Takes observation
    :return:
    """
    c = 299792428 # speed of light
    # print "=== Least Squares ==="
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
        l = np.matrix([P[i] - rho[i] + c*xyz[3] for i in xrange(len(sats))]).transpose()
        # form x-vector
        x_hat = ((AT*A).I * AT * l).flatten().getA()[0]
        x_hat[3] /= c
        # print "Q =",(AT*A).I.diagonal()
        # print "(x,y,z,cδt) =",x_hat
        # iterate
        xyz += x_hat
        delta = np.sqrt(sum(map(lambda k: k**2,x_hat[:3])+[x_hat[3]/c]))
        if delta < 10.:
            break
        # XYZs = np.array([navs[s][0].eph2pos(utc2gpst(o.date)) for s in sats])
        XYZs = np.array([navs[s][0].eph2pos(utc2gpst(o.date+dt.timedelta(seconds = x_hat[3]))) for s in sats])

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


    print lla_string(ecef_to_lat_lon_alt(least_squares(o, navigations)))

    # ecef_to_lat_lon_alt([3695041.,3695041.,3695041.])   # 45,45,0
    # ecef_to_lat_lon_alt([6400000*.707,6400000*.707,1])  # 0.45,0
    # ecef_to_lat_lon_alt([1,1,6400000])
    # ecef_to_lat_lon_alt([1,6400000,1])
    # ecef_to_lat_lon_alt([6400000,1,1])

    print sat_elev([6.38e6,0,0],[7.0e6,1e6,0])

    sat_positions, sat_names = [], []
    for s in navigations:
        n = navigations[s][0]
        xyz = n.eph2pos(n.date)
        sat_positions += [xyz]
        sat_names += [s]
        user_pos = least_squares(o, navigations)
        # print "User position:",ecef_to_lat_lon_alt(user_pos)
        print("Satellite's %s zenith angle: %.1f"%
              (s,sat_elev(user_pos,xyz)))

    user_pos = least_squares(o, navigations)[:3]
    print "User's position:",user_pos
    satellites(user_pos,sat_positions,sat_names)