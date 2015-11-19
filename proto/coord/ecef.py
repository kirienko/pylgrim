#! encoding: UTF8
import numpy as np
import gpstk as G
from numpy.linalg import norm
from numpy.matrixlib import matrix
from numpy.core.umath import sqrt, cos, sin, arcsin,isnan, degrees

# from proto.rover_position import xyz_string, lla_string


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


def ecef_to_lat_lon_alt1(R, deg = True):
    """
    Fukushima implementation of the Bowring algorithm (2006),
    see [4] --
    :param R: (X,Y,Z) -- coordinates in ECEF (numpy array)
    :return:  (φ,θ,h) -- lat [deg], lon [deg], alt in WGS84 (numpy array)
    """
    # WGS 84 constants
    a = 6378137.0       # Equatorial Radius [m]
    b = 6356752.314245  # Polar Radius [m]
    e = 0.08181919092890624 # e = sqrt(1-b²/a²)
    E = e**2
    e1 = sqrt(1 - e**2)            # e' = sqrt(1 - e²)
    if isinstance(R,list): R = np.array(R)
    p = sqrt(R[0]**2 + R[1]**2)             # (1) - sqrt(X² + Y²)
    az = abs(R[2])
    Z = e1*az/a
    P = p/a
    S, C = Z or 1, e1 * P or 1             # (C8) - zero approximation
    max_iter = 5
    for i in range(max_iter):
        A = sqrt(S**2 + C**2)
        Cn = P * A**3 - E * C**3
        Sn = Z * A**3 + E * S**3
        delta = abs(Sn/Cn - S/C)*C/S
        if isnan(delta):
            return ecef_to_lat_lon_alt1(R)
        if abs(delta) < 1e-10 or i == max_iter-1:
            break
        S, C = Sn, Cn
    # print "i =",i
    theta = np.math.atan2(R[1],R[0])
    Cc = e1 * Cn
    phi = np.sign(R[2])*np.math.atan2(Sn,Cc)
    h = (p*Cc + az*Sn - b*sqrt(Sn**2 + Cn**2))/sqrt(Cc**2 + Sn**2)
    if deg:
        out = np.array([np.degrees(phi), np.degrees(theta), h])
    else:
        out = np.array([phi, theta, h])
    # if filter(isnan, out):
    #     return ecef_to_lat_lon_alt1(R)
    return out


def ecef_to_lat_lon_alt(R, deg = True):
    """
    Fukushima implementation of the Bowring algorithm,
    see [3] -- equations (C7) - (C12)
    :param R: (X,Y,Z) -- coordinates in ECEF (numpy array)
    :return:  (φ,θ,h) -- lat [deg], lon [deg], alt in WGS84 (numpy array)
    """
    # WGS 84 constants
    a = 6378137.0       # Equatorial Radius [m]
    b = 6356752.314245  # Polar Radius [m]
    e = 0.08181919092890624 # e = sqrt(1-b²/a²)
    e1 = sqrt(1-e**2)            # e' = sqrt(1 - e²)
    c = a * e**2        # (6)
    if isinstance(R,list): R = np.array(R)
    p = sqrt(R[0]**2 + R[1]**2)             # (1) - sqrt(X² + Y²)
    T = R[2]/(e1*p) or 1.              # (C8) - zero approximation
    for i in range(10):
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
        h = sqrt(T1 - e**2)/e1 * (p - a/sqrt(T1))
    else:                                   # (C12)
        h = sqrt(T1 - e**2) * (R[2]/T - b/sqrt(T1))
    if deg:
        out = np.array([np.degrees(phi), np.degrees(theta), h])
    else:
        out = np.array([phi, theta, h])
    return out

def ecef_to_geodetic(R):
    """
    ECEF to Geodetic using gpstk. Degrees only
    """
    p = G.Position(*R)
    pos = np.array(G.Position.convertCartesianToGeodetic(p,6378137.0,0.006694380004260814))
    if pos[0] > 180: pos[0] = pos[0] - 360
    if pos[1] > 180: pos[1] = pos[1] - 360
    return pos

def sat_in_enu(R_u, R_sat):
    """
    Satellite coordinates in ENU at users position R_u
    :param R_u:   (X,Y,Z) -- user's coordinates in ECEF (numpy array)
    :param R_sat: (X,Y,Z) -- satellite's coordinates in ECEF (numpy array)
    :return: unit vector of the direction to the satellite in local `east, north, up` coords
    See [1] p. 522
    """
    # if isinstance(R_u,list):    R_u = np.array(R_u)
    if isinstance(R_sat,list):  R_sat = np.array(R_sat)
    phi, theta, h = ecef_to_lat_lon_alt(R_u, deg=False)
    # coordinate transformation matrix from ECEF to ENU:
    C_ecef_to_enu = matrix([[-sin(theta),cos(theta),0],
                            [-sin(phi)*cos(theta), -sin(phi)*sin(theta), cos(phi)],
                            [cos(phi)*cos(theta), cos(phi)*sin(theta),sin(phi)]])
    # coordinate origin shift vector from ECEF to local ENU:
    S_enu = matrix([[R_u[0] * sin(theta) - R_u[1] * cos(theta)],
                    [R_u[0] * sin(phi) * cos(theta) - R_u[1] * sin(phi) * sin(theta) - R_u[2] * cos(phi)],
                    [-R_u[0] * cos(phi) * cos(theta) - R_u[1] * cos(phi) * sin(theta) - R_u[2] * sin(theta)]])

    X_enu = C_ecef_to_enu * R_sat.reshape((3,1)) + S_enu
    # return normalized vector
    return (X_enu/norm(X_enu)).flatten().getA()[0]

def sat_elev(R_u, R_sat,deg = True):
    """
    Zenith angle of some certain sat. relative to the user's position
    :param R_u:   (X,Y,Z) -- user's coordinates in ECEF (numpy array)
    :param R_sat: (X,Y,Z) -- satellite's coordinates in ECEF (numpy array)
    :param deg: whether to return angle in degrees
    :return: elevation angle of the satellite
    """
    elev = arcsin(sat_in_enu(R_u,R_sat)[2])
    # print "Elev angle = %.1f deg" % elev
    if deg:
        return degrees(elev)
    else:
        return elev
