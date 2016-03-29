#! encoding: utf8

from unittest import TestCase, main
from proto.ils import qrmcp, reduction
from numpy.random import rand
from numpy import array, dot
from numpy.testing import assert_allclose, assert_array_equal

m, n = 5, 3
print "Initital size of matrix B: %d × %d" % (m, n)
# B = rand(m, n)
B = array([[-1.38483, 0.53704, 0.14925],  # ans:
           [1.05734, 0.61432, 0.94116],  # v1    v2
           [-0.33438, -0.13293, -0.60755],  # 1     1
           [0.26814, 0.41059, -0.52649],  # -2    -1
           [-0.66335, -1.42715, -0.97412]])  # 3     2
z_true = array([1, -2, 3]).reshape(3, 1)  # column vector


class TestQrmcp(TestCase):
    def test_qrmcp(self):
        y = dot(B, z_true) + 1e-3 * rand(m, 1)
        R_qrmcp, piv, y = qrmcp(B, y)
        R_true = [[-1.58217452, -1.20917892, -0.94591479],
                  [0, 1.19444506, -0.11444606],
                  [0, 0, 1.65879823]]
        assert_allclose(R_qrmcp, R_true, rtol=1e-7)


class TestReduction(TestCase):
    def test_reduction(self):
        y = dot(B, z_true) + 1e-3 * rand(m, 1)
        R_red, Z_red, y = reduction(B, y)
        R_true = [[ 1.25132918, -0.47161381, -0.39120129],
                  [ 0,           1.51025052,  0.86880047],
                  [-0,           0,          1.65879823]]
        Z_true = array([[0, 0, 1], [1, 0, 0], [-1, 1, 0]])
        assert_allclose(R_red, R_true, rtol=1e-5)
        assert_array_equal(Z_red.astype(int), Z_true)


if __name__ == "__main__":
     main()