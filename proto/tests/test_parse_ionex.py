from unittest import TestCase, main
from proto.helper.ionex import parse_ionex
import datetime as dt
import os


class TestParse_ionex(TestCase):
    def setUp(self):
        path = os.getcwd()
        if not os.path.basename(path) == 'data':
            os.chdir(os.path.join(path, 'proto', 'tests', 'data'))

    def test_parse_ionex(self):
        ionex = "igsg0010.16i"
        M = parse_ionex(ionex)
        dttimes = [dt.datetime(2016, 1, 1, x, 00) for x in range(24)]
        p = (60.1, 30.1)
        correct_vals = [0.913819144802, 0.896627150059, 0.782946952107, 0.643969126833,
                        0.55367705622, 0.553352311325, 0.617950565843, 0.906993007108,
                        1.10022271455, 1.33549389609, 1.66625306657, 1.60813672016,
                        1.4551818746, 1.28501554961, 0.919015063122, 0.570271520356,
                        0.420219894161, 0.389369129134, 0.388706649548, 0.291458543284,
                        0.243870426367, 0.268882278182, 0.342268129559, 0.374430863963]
        vals = [M(p, tt) for tt in dttimes]
        for i in range(len(vals)):
            TestCase.assertAlmostEqual(self, correct_vals[i], vals[i])


if __name__ == "__main__":
    main()
