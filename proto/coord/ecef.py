#! encoding: UTF8
import numpy as np

from numpy.core.umath import sqrt, cos, sin

from proto.rover_position import xyz_string, lla_string


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