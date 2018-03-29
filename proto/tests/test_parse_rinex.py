#! encoding: utf8

import datetime as dt
from unittest import TestCase, main
from proto.parse_rinex import parse_rinex, parse_sp3


class TestParse_rinex(TestCase):
    def test_parse_rinex_glo_nav(self):
        navigations = parse_rinex('../test_data/log_000.15g')
        g = navigations['R03']
        z1, z2 = sorted([g[0], g[1]], key=lambda x: x.date)
        t1, t2 = z1.eph[8], z2.eph[8]

        delta_t = dt.timedelta(seconds=t2-t1)
        TestCase.assertAlmostEqual(self, delta_t.total_seconds(), 5469.878906, places=10)

    # TODO: implement the rest:
    # def test_parse_rinex_gps_nav(self):
    # def test_parse_rinex_obs(self):

    def test_parse_rinex_sp3(self):
        precise = parse_sp3('../test_data/igs11484.sp3')
        TestCase.assertEqual(self, len(precise), 28)
        TestCase.assertEqual(self, precise['G01'][0].date.std.year, 2002)


if __name__ == '__main__':
    main()
