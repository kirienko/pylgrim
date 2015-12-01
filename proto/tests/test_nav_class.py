from numpy import array, sqrt
from unittest import TestCase
from numpy.testing import assert_allclose
from proto.nav_data import *
from proto.parse_rinex import parse_rinex
from datetime import timedelta


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
        nav_file = '../../test_data/epgga2.010'
        navigations = parse_rinex(nav_file)

        self.t = timedelta(seconds=346500)

        self.n = navigations['G31'][0]
        self.right_pos = array([11660.379642, 11313.211213, -21326.822815])

    def test_eph2pos(self):
        XYZ = self.n.eph2pos(self.n.epoch + self.t) / 1000
        assert_allclose(XYZ, self.right_pos, rtol=1e-6)
        dist = sqrt(sum(map(lambda x, y: (x - y) ** 2, XYZ, self.right_pos)))
        assert_allclose(dist, 0., atol=1e-2)
        print "%.2f m <-- distance between exact and calculated positions" % (dist * 1000)


class TestNavGLO(TestCase):
    def setUp(self):
        nav_file = '../../test_data/test.g'
        navigations = parse_rinex(nav_file)
        self.navs = navigations['R19']
        # for j, n in enumerate(self.navs):
        #     print j, n.date
        # print self.navs[1].date - self.navs[0].date

    @staticmethod
    def middle(x):
        """
        :param x: array of length N
        :return: array of length N-1 of the middles between elemes of x
        """
        return map(lambda x: x[0].date +
                             timedelta(seconds=float((x[1].date - x[0].date).total_seconds()) / 4),
                   zip(x[:-1], x[1:]))

    def test_eph2pos(self):
        dt = self.middle(self.navs)
        # print "dt =", dt
        d1 = self.navs[2].eph2pos(dt[1])
        d2 = self.navs[1].eph2pos(dt[1])
        # print "d1 =", d1
        # print "d2 =", d2
        assert_allclose(d1, d2, rtol=1e-6)
