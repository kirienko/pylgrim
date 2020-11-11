import os
from numpy import array, sqrt
from random import choice
from unittest import TestCase, skip
from numpy.testing import assert_allclose
from proto.nav_data import *
from proto.nav_data import def_leap
from proto.parse_rinex import parse_rinex
from datetime import datetime


def distance(R1, R2):
    """
    Calculates euclidean distance (along the straight line)
    :param R1: vector in ECEF
    :param R2: vector in ECEF
    :return: Euclidean distance between R1 and R2 [in meters]
    """
    return sqrt(sum(map(lambda x, y: (x - y) ** 2, R1, R2)))


class TestNavGPS(TestCase):
    """
    Source: http://web.ics.purdue.edu/~ecalais/teaching/gps_geodesy/lab_4.pdf
    """

    def setUp(self):
        path = os.getcwd()
        if not os.path.basename(path) == 'data':
            os.chdir(os.path.join(path, 'proto', 'tests', 'data'))
        nav_file = 'epgga2.010'
        navigations = parse_rinex(nav_file)

        self.t = 346500     # seconds from epoch
        self.n = navigations['G31'][0]
        self.right_pos = array([11660.379642, 11313.211213, -21326.822815])

    def test_eph2pos(self):
        XYZ = self.n.eph2pos(self.n.epoch + self.t) / 1000
        dist = sqrt(sum(map(lambda x, y: (x - y) ** 2, XYZ, self.right_pos)))
        print("%.2f m <-- distance between exact and calculated positions" % (dist * 1000))
        assert_allclose(XYZ, self.right_pos, rtol=1e-6)
        assert_allclose(dist, 0., atol=1e-2)

    @skip("Not implemented")
    def test_utc2gps(self):
        self.fail("utc2gps test is not implemented")


class TestNavGLO(TestCase):
    def setUp(self):
        path = os.getcwd()
        if not os.path.basename(path) == 'data':
            os.chdir(os.path.join(path, 'proto', 'tests', 'data'))
        nav_file = 'log_000.15g'
        navigations = parse_rinex(nav_file)
        sats = [k for k, v in navigations.items() if len(v) > 1]   # list of possible satellites
        sat = choice(sats)  # choose one satellite randomly
        print("PRN: " + sat)
        self.navs = navigations[sat]
        # self.navs = navigations['R19']
        # for j, n in enumerate(self.navs):
        #     print j, n.date
        # print self.navs[1].date - self.navs[0].date

    @staticmethod
    def middle(x):
        """
        :param x: array of length N
        :return: array of length N-1 of the middles between elements of x
        """
        return list(map(lambda z: z[0].date + ((z[1].date - z[0].date) / 2),
                   zip(x[:-1], x[1:])))

    def test_eph2pos(self):
        dt = self.middle(self.navs)
        print("dt =", dt)
        num = 0
        d1 = self.navs[num].eph2pos(dt[num])
        d2 = self.navs[num+1].eph2pos(dt[num])
        print("at {} d1 = {}".format(dt[num], d1, self.navs[num].date))
        print("at {} d2 = {}".format(dt[num], d2, self.navs[num+1].date))
        assert_allclose(d1, d2, rtol=1e-2)  # FIXME: usually it fails


class TestDef_leap(TestCase):
    def test_def_leap(self):
        self.assertEqual(def_leap(datetime(1980, 11, 22, 0, 0, 0)), 0)
        self.assertEqual(def_leap(datetime(1981, 11, 22, 0, 0, 0)), 1)
        self.assertEqual(def_leap(datetime(1982, 11, 22, 0, 0, 0)), 2)
        self.assertEqual(def_leap(datetime(1983, 11, 22, 0, 0, 0)), 3)
        self.assertEqual(def_leap(datetime(1985, 11, 22, 0, 0, 0)), 4)
        self.assertEqual(def_leap(datetime(1988, 11, 22, 0, 0, 0)), 5)
        self.assertEqual(def_leap(datetime(1990, 11, 22, 0, 0, 0)), 6)
        self.assertEqual(def_leap(datetime(1991, 11, 22, 0, 0, 0)), 7)
        self.assertEqual(def_leap(datetime(1992, 11, 22, 0, 0, 0)), 8)
        self.assertEqual(def_leap(datetime(1993, 11, 22, 0, 0, 0)), 9)
        self.assertEqual(def_leap(datetime(1994, 11, 22, 0, 0, 0)), 10)
        self.assertEqual(def_leap(datetime(1996, 11, 22, 0, 0, 0)), 11)
        self.assertEqual(def_leap(datetime(1998, 11, 22, 0, 0, 0)), 12)
        self.assertEqual(def_leap(datetime(1999, 11, 22, 0, 0, 0)), 13)
        self.assertEqual(def_leap(datetime(2006, 11, 22, 0, 0, 0)), 14)
        self.assertEqual(def_leap(datetime(2009, 11, 22, 0, 0, 0)), 15)
        self.assertEqual(def_leap(datetime(2012, 11, 22, 0, 0, 0)), 16)
        self.assertEqual(def_leap(datetime(2015, 11, 22, 0, 0, 0)), 17)