#! encoding: utf8
from math import cos, sin, asin, tan, pi, sqrt

D2R = pi / 180.0


def ipp(pos, azel, r_earth, hion):
    """
    Compute ionospheric pierce point (ipp) position and slant factor
    :param pos: receiver position {lat,lon,h} (rad,m)
    :param azel: azimuth/elevation angle {az,el} (rad)
    :param r_earth: earth radius (km)
    :param hion: altitude of ionosphere (km)
    :return: pierce point position {lat,lon,h} (rad,m), slant factor
    """
    rp = r_earth / (r_earth + hion) * cos(azel[1])
    ap = pi / 2.0 - azel[1] - asin(rp)
    sinap = sin(ap)
    tanap = tan(ap)
    cosaz = cos(azel[0])
    posp = []
    posp += [asin(sin(pos[0]) * cos(ap) + cos(pos[0]) * sinap * cosaz)]

    if (pos[0] > 70.0 * D2R and tanap * cosaz > tan(pi / 2.0 - pos[0])) or \
            (pos[0] < -70.0 * D2R and -tanap * cosaz > tan(pi / 2.0 + pos[0])):
        posp += [pos[1] + pi - asin(sinap * sin(azel[0]) / cos(posp[0]))]

    else:
        posp += [pos[1] + asin(sinap * sin(azel[0]) / cos(posp[0]))]

    return posp, 1.0 / sqrt(1.0 - rp * rp)
