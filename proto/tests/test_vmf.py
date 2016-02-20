#! encoding: utf8
from datetime import datetime
from numpy import deg2rad
from numpy.testing import assert_allclose
from unittest import TestCase

from proto.vmf import find_VMF_coeffs, vmf


class TestVmf(TestCase):
    def setUp(self):
        self.boehm_time = datetime(2008, 03, 19, 0, 0)
        self.boehm_pos_deg = (46, 15, 1000)
        self.boehm_pos_rad = (deg2rad(46), deg2rad(15), 1000)
        self.ahf = "../../test_data/ah08079.h00"   # ah = 0.00121328
        self.awf = "../../test_data/aw08079.h00"   # aw = 0.00043331
        self.boehm_coeffs = (10.1804, 10.9260)
        self.cfs = find_VMF_coeffs(self.ahf, self.awf, self.boehm_pos_deg)

    def test_find_VMF_coeffs(self):
        self.assertEquals(self.cfs, (0.00121328, 0.00043331))

    def test_vmf(self):
        ori = vmf(pos=self.boehm_pos_rad, time=self.boehm_time, coeffs=self.cfs, zd=deg2rad(85))
        assert_allclose(ori, self.boehm_coeffs, rtol=1e-5, err_msg="Wrong tropodelays")
