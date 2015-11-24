from unittest import TestCase
from proto.nav_data import def_leap
from datetime import datetime


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
