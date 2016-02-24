#!/usr/bin/python
# ! encoding: UTF8
from math import cos, exp, pi, radians
from datetime import datetime

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

if __name__ == "__main__":
    pos = [pi / 3, pi / 6, 100.]
    azel = 45 * pi / 180
    t = datetime(2016, 1, 1, 0, 0)
    from vmf import vmf

    vmf_coeffs = (0.00121328, 0.00043331)
    for h in range(10):
        saa_ = tropmodel(pos, azel, 0.1 * h)
        print "Relative humidity: %2d %% \t delay: %.3f meters" % (h * 10, saa_)
    for h in range(90, 0, -10):
        saa_ = tropmodel(pos, radians(h), 50)
        vmf_ = sum(vmf(pos, t, radians(h),  vmf_coeffs))
        print "Elevation angle: %2d° \t Delays [m]: Saastamoinen: %.3f, VMF1: %.3f" % (h, saa_, vmf_)
