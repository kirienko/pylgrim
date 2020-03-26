#! encoding: utf8
from __future__ import absolute_import
import os
from unittest import TestCase
from numpy.testing import assert_allclose
from proto.parse_rinex import parse_rinex, parse_sp3
from proto.least_squares import least_squares
from proto.helper.vmf import find_VMF_coeffs
from proto.coord.ecef import ecef_to_lat_lon_alt as ecef2lla


class TestLeast_squares(TestCase):
    """
    Test case from here: http://web.ics.purdue.edu/~ecalais/teaching/gps_geodesy/lab_6.pdf  # noqa
    """
    def setUp(self):
        td_path = os.path.join(os.getcwd(), 'proto', 'tests', 'data')
        os.chdir(td_path)
        self.apriori_coords = [4433070.0, 362070.0, 4556010.0]
        self.observations = parse_rinex('sjdv0100.02o')
        self.sp3s = parse_sp3('igs11484.sp3')
        self.o = self.observations[29]
        print(self.o.date, "<-- okay, that's it")
        self.vmf_coeffs = find_VMF_coeffs('ah02010.h00', 'aw02010.h00',
                                          ecef2lla(self.apriori_coords))
        # self.vmf_coeffs = (0.00122794, 0.00046177)
        print("VMF coefficients:", self.vmf_coeffs)

    def test_least_squares(self):
        print("PRN number:", self.o.PRN_number)
        # make sure that GPS satellites are used (not GLONASS)
        pos = least_squares(self.o, self.sp3s, init_pos=self.apriori_coords, vmf_coeffs=self.vmf_coeffs)
        approx_ans = [4433432, 362722, 4556149]
        print(pos, approx_ans)
        assert_allclose(pos, approx_ans, rtol=1e-5)
