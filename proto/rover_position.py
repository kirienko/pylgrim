#!/usr/bin/python
#! encoding: UTF8

import datetime as dt
import numpy as np

from coord.ecef import ecef_to_lat_lon_alt, sat_elev
from parse_rinex import parse_rinex
from visualization.ellipsoid import satellites
from visualization.map import on_map
from delays import tropmodel

__author__ = 'kirienko'

"""
1. Parse Nav file --> ephemerids --> sat coords (ECEF)
2. Parse Obs file --> pseudoranges to rover
3. Calculate rover's coords in ECEF
"""

# Pretty printing
def lla_string(R):
    return "φ = %.3f°, θ = %.3f°, h = %d m" % (R[0],R[1],int(R[2]))
def xyz_string(R):
    return "(%d, %d, %d) [km]" % tuple(map(int,R / 1000))

def nav_nearest_in_time(t,nav_array):
    '''
    From array of NavGPS objects returns the one
        which ephemeris are the closets in time to ``t``
    :param t: UTC time
    :param nav_array:
    :return:
    '''
    diff_array = [abs((n.date - n.utc2gps(t)).total_seconds()) for n in nav_array]
    return nav_array[diff_array.index(min(diff_array))]

def least_squares(obs, navs, init_pos = ''):
    """
    x = (A^TA)^{-1}A^T l
    Takes an observation ``obs`` and all the data ``nav`` from navigation file.
    If we have a-priori information about rover's position,
        then we can filter low satellites and use troposperic correction
    :return: rover's position in ecef [m]
    """
    c = 299792428   # speed of light
    elev_mask = 10  # satellite elevation mask
    now = obs.date
    # Find all possible satellites N
    sats = []
    # sats = {}
    for i,r in enumerate(obs.PRN_number):
        if obs.obs_data['C1'][i] and obs.obs_data['P2'][i] and 'G' in r:
            nnt = nav_nearest_in_time(now,navs[r])
            if len(init_pos):
                sat_coord = nnt.eph2pos(now)
                if sat_elev(init_pos,sat_coord) < elev_mask:
                    print "Satellite %s excluded" % r
                    continue
            # sats += [r]
            sats += [(r,nnt)]
            # sats.update({r:nnt})
    # Form matrix if N >= 4:
    if len(sats) > 3:
        # observed [iono-free] pseudoranges
        P = np.array([obs.ionofree_pseudorange(s[0]) for s in sats])
        # print "P =",P
        # get XYZ-coords of satellites
        # XYZs = np.array([nav_nearest_in_time(now,navs[s]).eph2pos(now) for s in sats])
        XYZs = np.array([s[1].eph2pos(now) for s in sats])
        # print "XYZs =",XYZs
    elif len(sats) <= 3 and len(init_pos):     # FIXME: rewise this logic
        print "\n\tWarning: too few satellites:", len(sats)
        return None
    else:
        print "\n\tWarning: bad measurement!"
        return None
    # if err == {}: err = {s[0]:0. for s in sats}
    xyzt = [1e-10,1e-10,1e-10,0.] # initial point
    if len(init_pos):
        xyzt = init_pos + [0.]
    for itr in range(10):
        # print "\t iter =", itr,
        # geometrical ranges
        lla = ecef_to_lat_lon_alt(xyzt)
        rho = np.array([np.sqrt(sum([(x - xyzt[i])**2 for i,x in enumerate(XYZs[j])])) for j in xrange(len(sats))])
        # from A-matrix
        A = np.matrix([np.append((xyzt[:3] - XYZs[i])/rho[i],1) for i in xrange(len(sats))])
        AT = A.transpose()
        # form l-vector (sometimes `l` is denoted as `b`)
        l = np.matrix([P[i] - rho[i] + c*s[1].time_offset(now) - tropmodel(lla,sat_elev(xyzt[:3], XYZs[i]))
                       for i,s in enumerate(sats)]).transpose()
        # form x-vector
        x_hat_matrix = ((AT*A).I * AT * l)
        x_hat = x_hat_matrix.flatten().getA()[0]
        x_hat[3] /= c
        # print "(x,y,z,cδt) =",x_hat
        xyzt += x_hat
        # print lla_string(ecef_to_lat_lon_alt(xyzt)),"%.4f"%xyzt[3]
        delta = np.sqrt(sum(map(lambda k: k**2,x_hat[:3])))
        if delta < 10.:
            break
        # now += dt.timedelta(seconds=x_hat[3])
        XYZs = np.array([s[1].eph2pos(now) for s in sats])

    if len(init_pos):
        # Q = (AT*A).I.diagonal().tolist()[0]
        # print "Horizontal DOP: %.2f m" % np.sqrt(Q[0]**2 + Q[1]**2)
        # print l - A*x_hat_matrix
        return xyzt
    else:
        # errors = {s[0]:(l - A*x_hat_matrix).tolist()[i][0] for i,s in enumerate(sats)}
        # print errors
        # print "try with initial position",xyzt,
        return least_squares(obs, navs, xyzt)

if __name__ == "__main__":

    nav_file = '../test_data/test.n'
    obs_file = '../test_data/test.o'

    # Process Nav file:
    # ``navigations`` is a dict with
    #       keys:   GNSS identificators, i.e. 'G16', 'G04', ...
    #       values: Nav observation objects
    #   Note: One satellite may have several nav objects (for several times,
    #       e.g. data on 14:00 and on 16:00)
    navigations = parse_rinex(nav_file)

    # Process Obs file
    observations = parse_rinex(obs_file)
    o = observations[240]
    print o.sat_types



    '''
    sat_positions, sat_names = [], []
    user_pos = least_squares(o, navigations)
    for s in navigations:
        n = navigations[s][0]
        xyz = n.eph2pos(n.date)
        sat_positions += [xyz]
        sat_names += [s]
        # print "User position:",ecef_to_lat_lon_alt(user_pos)
        print("Satellite's %s zenith angle: %.1f"%
              (s,sat_elev(user_pos,xyz)))
    '''
    # satellites(user_pos,sat_positions,sat_names)
    user_pos = []
    for num_o in range(190,250,10):
        print
        user_pos += [least_squares(observations[num_o], navigations)]
    user_pos = [up[:3] for up in user_pos if up is not None]
    delta = lambda x,y: (x - y)**2
    dd =  map(delta,user_pos[1:],user_pos[:-1])
    dd = [int(np.sqrt(sum(up))) for up in dd]
    print dd
    print "User's position:\n",'\n'.join(map(lambda x: lla_string(ecef_to_lat_lon_alt(x)),user_pos))
    # on_map(map(ecef_to_lat_lon_alt,user_pos))
