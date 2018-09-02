#! encoding: utf8

from unittest import TestCase, main
from proto.ils import qrmcp, reduction, ils
from numpy.random import rand, randint
from numpy import array, dot
from numpy.testing import assert_allclose, assert_array_equal


m, n = 5, 3
print("Initital size of matrix B: %d × %d" % (m, n))
# B = rand(m, n)
B = array([[-1.38483, 0.53704, 0.14925],        # ans:
           [1.05734, 0.61432, 0.94116],         # v1    v2
           [-0.33438, -0.13293, -0.60755],      # 1     1
           [0.26814, 0.41059, -0.52649],        # -2    -1
           [-0.66335, -1.42715, -0.97412]])     # 3     2
z_true = array([1, -2, 3]).reshape(3, 1)        # column vector


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
        R_true = [[ 1.95477311, -0.77764188, -1.12307502],
                  [ 0,          -1.89700789, -1.29743768],
                  [-0,           0,          -1.86127017]]
        Z_true = array([[1, 0, 0], [-1, 1, 0], [0, 0, 1]])
        assert_allclose(R_red, R_true, rtol=1e-5)
        assert_array_equal(Z_red.astype(int), Z_true)


class TestIls(TestCase):
    def test_ils(self):
        for j in range(10):
            zt = randint(low=-10, high=10, size=(3, 1))
            y = dot(B, zt) + 1e-3 * rand(m, 1)
            Z = ils(B.copy(), y.copy())
            assert_array_equal(Z, zt)


if __name__ == "__main__":
     main()


# well, I don't think we need it anymore
def pprint_two_matrices(M1, M2, dtype='float', name1='', name2=''):
    """
    Pretty print two matrices (of the same size) in parallel
    """
    if len(M1) != len(M2):
        if len(M1.reshape(1, max(M1.shape))) != len(M2.reshape(1, max(M2.shape))):
            print("Matrices are of different length: %d and %d" % (len(M1), len(M2)))
            return
        else:
            _M1 = M1.reshape(1, max(M1.shape)).astype(dtype)
            _M2 = M2.reshape(1, max(M2.shape)).astype(dtype)
    else:
        _M1, _M2 = M1.astype(dtype), M2.astype(dtype)
    len_M = len(_M1)
    len_str = lambda x: len(str(x))
    max_elem_len = min(8, max(map(len_str, M1.flatten()) + map(len_str, M2.flatten())))
    width = (len(_M1[0])) * (max_elem_len + 3)
    eq = "=" if (_M1 == _M2).all() else "≠"
    print()
    if name1 or name2:
        print("{:^{}}".format(name1, width) + eq + "{:^{}}".format(name2, width))
    for j in xrange(len_M):
        q1, q2 = map(str, _M1[j]), map(str, _M2[j])
        str1 = " ".join(map(lambda x: '{:>{}.{}}'.format(x, max_elem_len + 2, max_elem_len), q1))
        str2 = " ".join(map(lambda x: '{:>{}.{}}'.format(x, max_elem_len + 2, max_elem_len), q2))
        print(str1 + ' | ' + str2)