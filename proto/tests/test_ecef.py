from unittest import TestCase, main
from proto.coord.ecef import *
from numpy.testing import assert_allclose
from numpy import array


class TestEcef_to_lat_lon_alt(TestCase):
    def setUp(self):
        self.inp = array([[-576793.17, -5376363.47, 3372298.51], [2297292.91, 1016894.94, -5843939.62]])
        self.ans = array([[32.12345, -96.12345, 500.0], [-66.87654, 23.87654, 1000.0]])

    def test_ecef_to_lat_lon_alt(self):
        lla1, lla2 = map(ecef_to_lat_lon_alt, self.inp)
        assert_allclose(lla1, self.ans[0], rtol=1e-5, err_msg="Wrong ECEF to LatLon conversion")
        assert_allclose(lla2, self.ans[1], rtol=1e-5, err_msg="Wrong ECEF to LatLon conversion")

    def test_ecef_to_lat_lon_alt1(self):
        lla1, lla2 = map(ecef_to_lat_lon_alt1, self.inp)
        assert_allclose(lla1, self.ans[0], rtol=1e-5, err_msg="Wrong ECEF to LatLon conversion")
        assert_allclose(lla2, self.ans[1], rtol=1e-5, err_msg="Wrong ECEF to LatLon conversion")


if __name__ == "__main__":
    main()
