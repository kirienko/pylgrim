#!/usr/bin/python
# ! encoding: UTF8
from datetime import datetime
from math import sin, cos, exp, pi, radians, floor

from gtime import GTime

__author__ = 'kirienko'


## based on RTKlib

def saast(pos, el, humi=0.75, temp0=15.0):
    """
    Function from RTKlib: https://github.com/tomojitakasu/RTKLIB/blob/master/src/rtkcmn.c#L3362-3362
        with no changes
    :param time:    time
    :param pos:     receiver position {lat,lon,h} (rad,m)
    :param el:    azimuth/elevation angle {az,el} (rad) -- we do not use az
    :param humi:    relative humidity
    :param temp0:   temperature (Celsius)
    :return:        tropospheric delay (m)
    """
    # temp0 = 15.0  # temparature at sea level
    # double hgt,pres,temp,e,z,trph,trpw;

    # if (pos[2]<-100.0||1E4<pos[2]||el[1]<=0) return 0.0;
    if pos[2] < -100.0 or 1E4 < pos[2] or el <= 0:
        return 0.0

    # /* standard atmosphere */
    # hgt=pos[2]<0.0?0.0:pos[2];
    hgt = 0.0 if pos[2] < 0.0 else pos[2]

    pres = 1013.25 * pow(1.0 - 2.2557E-5 * hgt, 5.2568)
    temp = temp0 - 6.5E-3 * hgt + 273.16
    e = 6.108 * humi * exp((17.15 * temp - 4684.0) / (temp - 38.45))

    # /* saastamoninen model */
    z = pi / 2.0 - el
    trph = 0.0022768 * pres / (1.0 - 0.00266 * cos(2.0 * pos[0]) - 0.00028 * hgt / 1E3) / cos(z)
    trpw = 0.002277 * (1255.0 / temp + 0.05) * e / cos(z)
    return trph + trpw


def tropmodel(pos, el, time='', coeffs=(), humi=0.75, temp0=15.0):
    if len(coeffs) > 1:
        return vmf(pos, time, el, coeffs)
    else:
        return saast(pos, el, humi, temp0)



def klobuchar(pos, az, el, time, iono_coeffs):
    """
    Details are taken from [5]: IS-GPS-200H, Fig. 20-4
    Note: result is referred to the GPS L₁ frequency;
    if the user is operating on the GPS L₂ frequency, the correction term must
    be multiplied by γ = f₂²/f₁¹ = 0.6071850227694382
    :param pos: [lat, lon, alt] in radians and meters
    """

    if pos[2] < -1E3 or el < 0:
        return 0.0
    if len(iono_coeffs) < 8:
        return None

    # earth centered angle (semi-circle) 
    psi = 0.0137 / (el / pi + 0.11) - 0.022

    # subionospheric latitude/longitude (semi-circle)
    phi = pos[0] / pi + psi * cos(az)
    if phi > 0.416:
        phi = 0.416
    elif phi < -0.416:
        phi = -0.416
    lam = pos[1] / pi + psi * sin(az) / cos(phi * pi)

    # geomagnetic latitude (semi-circle) */
    phi += 0.064 * cos((lam - 1.617) * pi)

    # local time (s)
    tt = 43200.0 * lam + time2gpst(time)
    tt -= floor(tt / 86400.0) * 86400.0     # 0<=tt<86400

    # slant factor
    f = 1.0 + 16.0 * pow(0.53 - el / pi, 3.0)

    # ionospheric delay 
    amp = iono_coeffs[0] + phi * (iono_coeffs[1] + phi * (iono_coeffs[2] + phi * iono_coeffs[3]))
    per = iono_coeffs[4] + phi * (iono_coeffs[5] + phi * (iono_coeffs[6] + phi * iono_coeffs[7]))
    if amp < 0.0:
        amp = 0.
    if per < 72000.0:
        per = 72000.0
    x = 2.0 * pi * (tt - 50400.0) / per

    if abs(x) < 1.57:
        return 2.99792458E8 * f * (5E-9 + amp * (1.0 + x * x * (-0.5 + x * x / 24.0)))
    else:
        return 2.99792458E8 * f * 5E-9


def time2gpst(t):
    """
    convert GTime to time of week in gps time
    :param t: GTime instance
    :return: time of week in gps time (s)
    """
    t0 = GTime(1980, 1, 6, 0, 0, 0)  # GPS time reference
    sec = t - t0  # number of seconds since t0, type: float
    gps_week = int(sec / (86400 * 7))

    return (sec - gps_week * 86400 * 7) + (t.sec - int(t.sec))


if __name__ == "__main__":
    pos = [pi / 3, pi / 6, 100.]
    azel = 45 * pi / 180
    t = datetime(2016, 1, 1, 0, 0)
    from proto.helper.vmf import vmf

    vmf_coeffs = (0.00121328, 0.00043331)
    for h in range(10):
        saa_ = tropmodel(pos, azel, 0.1 * h)
        print "Relative humidity: %2d %% \t delay: %.3f meters" % (h * 10, saa_)
    for h in range(90, 0, -10):
        saa_ = tropmodel(pos, radians(h), 50)
        vmf_ = sum(vmf(pos, t, radians(h),  vmf_coeffs))
        print "Elevation angle: %2d° \t Delays [m]: Saastamoinen: %.3f, VMF1: %.3f" % (h, saa_, vmf_)
