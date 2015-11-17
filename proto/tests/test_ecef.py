from unittest import TestCase, main
from proto.coord.ecef import ecef_to_lat_lon_alt


class TestEcef_to_lat_lon_alt(TestCase):
    def test_ecef_to_lat_lon_alt(self):
        lla1 = ecef_to_lat_lon_alt([-576793.17, -5376363.47, 3372298.51])
        ans1 = (32.12345, -96.12345, 500.0)
        for i in xrange(3):
            self.assertAlmostEqual(lla1[i], ans1[i])

        lla2 = ecef_to_lat_lon_alt([2297292.91, 1016894.94, -5843939.62])
        ans2 = (-66.87654, 23.87654, 1000.0)
        for i in xrange(3):
            self.assertAlmostEqual(lla2[i], ans2[i])

if __name__ == "__main__":
    main()