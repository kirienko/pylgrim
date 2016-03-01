from unittest import TestCase
from proto.gtime import GTime


class TestGTime(TestCase):
    def setUp(self):
        self.b = GTime(2000, 1, 1, 0, 2, 3.45)
        self.z = GTime(2000, 1, 1, 0, 0, 10.45)

    def test__init__(self):
        self.assertEqual(self.b.std.timetuple()[:5], (2000, 1, 1, 0, 2))
        self.assertEqual(self.z.std.timetuple()[:5], (2000, 1, 1, 0, 0))
        self.assertEqual(self.b.sec, 3.45)
        self.assertEqual(self.z.sec, 10.45)

    def test__sub__(self):
        self.assertAlmostEqual(self.b - GTime(2000, 1, 1, 0, 2, 1.05), 2.4, 10)
        self.assertAlmostEqual(self.b - GTime(2000, 1, 1, 0, 0, 0), 123.45, 10)
        self.assertAlmostEqual(self.z - self.b, -113., 10)

    def test__add__(self):
        t = list()
        t.append(self.z + 10)                           # 2000-01-01 00:00:20.450000000
        t.append(self.z + 50.55)                        # 2000-01-01 00:01:01.000000000
        t.append(self.z + 3650.55)                      # 2000-01-01 01:01:01.000000000
        t.append(self.z + 3650.55 + 60 * 60 * 24 * 2)   # 2000-01-03 01:01:01.000000000
        self.assertSequenceEqual(map(lambda x: x.std.year, t), [2000] * len(t))
        self.assertSequenceEqual(map(lambda x: x.std.month, t), [1] * len(t))
        self.assertSequenceEqual(map(lambda x: x.std.day, t), [1, 1, 1, 3])
        self.assertSequenceEqual(map(lambda x: x.std.hour, t), [0, 0, 1, 1])
        self.assertSequenceEqual(map(lambda x: x.std.minute, t), [0, 1, 1, 1])
        self.assertIsInstance(t[0].sec, float)
