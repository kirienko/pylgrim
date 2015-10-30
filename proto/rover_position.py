#!/usr/bin/python
#! encoding: UTF8

import datetime as dt
from parse_rinex import parse_rinex, gps2utc, LEAP, A0, A1
import numpy as np

__author__ = 'kirienko'

"""
1. Parse Nav file --> ephemerids --> sat coords (ECEF)
2. Parse Obs file --> pseudoranges to rover
3. Calculate rover's coords in ECEF
"""

def dummy_rover_pos(rho,R):
    """
    Naïve approach explained in [1], 2.4 (p. 48):
    calculate rover's position using measurement of 4 pseudoranges to 4 satellites.
    :param: r = (r1,r2,r3,r4) -- pseudoranges to sat-s 1,2,3,4
    :param: R = (R1,R2,R3,R4) -- coordinates  of sat-s 1,2,3,4 at the moment of measurement
            Ri = (Xi,Yi,Zi) -- coordinates of i-th satellite in ECEF
    :return: X = M⁻¹ × Y
    """
    sq = lambda x: np.square(x).sum()
    Y = np.array([rho[i]-sq(R[i])-ear for i in range(4)])
    M = np.matrix([np.append(-2*R[i],1) for i in range(4)])
    return M.I*Y

if __name__ == "__main__":

    nav_file = '../test_data/test.n'
    obs_file = '../test_data/test.o'

    # Process Nav file:
    # ``navigations`` is a dict with
    #       keys:   GNSS identificators, i.e. 'G16', 'G04', ...
    #       values: Nav observation objects
    navigations = parse_rinex(nav_file)

    # Process Obs file
    # observations = parse_rinex(obs_file)


    print "\nSatellie side:"
    satellites = navigations.keys()
    # print satellites
    print "Assume %s:" % satellites[0]

    o1 = navigations[satellites[0]][0]
    o2 = navigations[satellites[0]][1]
    delta_t = dt.timedelta(hours=1)
    print o2.eph2pos(o2.date-delta_t)
    print o1.eph2pos(o1.date+delta_t)