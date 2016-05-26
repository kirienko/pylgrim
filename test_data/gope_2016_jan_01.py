#! encoding: utf8

"""
Here we test the absolute accuracy of our software in static measurements.
GOPE = geodetic observatory Pecný: http://www.pecny.cz/gop/index.php/gnss/observations/gope
"""

from proto.coord.ecef import ecef_to_lat_lon_alt as ecef2lla
from proto.helper.vmf import find_VMF_coeffs
from proto.parse_rinex import parse_rinex, parse_sp3
from proto.least_squares import least_squares, distance

if __name__ == '__main__':
    home = './'
    sp3_file = home + 'igs18775.sp3'    # <-- GPS
    # sp3_file = home + 'igl18775.sp3'    # <-- GLONASS
    obs_file = home + 'gope0010.16o.GPS.filtered'
    # obs_file = home + 'gope0010.16o.GLO.filtered'
    navs = parse_sp3(sp3_file)
    obs  = parse_rinex(obs_file)

    zero_epoch = obs[2]
    print zero_epoch.date
    a_priori_ecef = [3979316.4389, 1050312.2534, 4857066.9036]  # GOPE coordinates in meters
    coords = ecef2lla(a_priori_ecef)
    print "A priori coordinates:", coords
    print "A priori ECEF coords:", a_priori_ecef
    vmf = find_VMF_coeffs(home+'VMF_ah16001.h00', home+'VMF_aw16001.h00', coords)
    print "VMF1 coefficients: ah = %.7f, aw = %.7f" % vmf
    for o in obs:
        final = least_squares(obs=o, navs=navs, init_pos=a_priori_ecef, vmf_coeffs=vmf)
        print str(o.date)[11:19], "Accuracy: %.1f m" % distance(a_priori_ecef, final)