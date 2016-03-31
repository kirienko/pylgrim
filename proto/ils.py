#! encoding: utf8

"""
Integer least squares solver,
MILES reimplemented in python.
See:    Chang, Xiao-Wen, and Tianyang Zhou.
        "MILES: MATLAB package for solving Mixed Integer LEast Squares problems."
        GPS Solutions 11.4 (2007): 289-294.
        DOI: 10.1007/s10291-007-0063-y
"""
from __future__ import division
from numpy import arange, array, dot, inf, ones, sqrt, triu, zeros
from numpy.linalg import lstsq, norm, solve, qr
from numpy.linalg import matrix_rank as rank
# for tests:
from oct2py import octave
from os.path import expanduser as eu
octave.addpath(eu('~')+'/code/pylgrim/ils')

def qrmcp(B, y):
    """
    QR factorization of B.
    Input arguments:
    B - m by n real matrix to be factorized (np.array)
    y - m-dimensional real vector to be transformed to Q'y

 Output arguments:
    R   - n by n real upper triangular matrix
    piv - n-vector storing the information of the permutation matrix P
    y   - m-vector transformed from the input y by Q, i.e., y := Q'*y
    """
    # Check input arguments
    m, n = B.shape
    if m < n:
        raise ValueError("Matrix to be factorized is column-rank deficient!")

    m2, n2 = y.shape
    if m != m2 and n2 != 1:
        raise ValueError("Input arguments have a dimension error!")

    # Initialization
    colnormB = zeros((2, n))
    piv = arange(n)

    # Compute the 2-norm squared of each column of B
    for j in xrange(n):
        colnormB[0][j] = norm(B[:, j]) ** 2

    for k in xrange(n):
        # Find the column with minimum 2-norm in B(k+1:m,k+1:n)
        tmp = colnormB[0, k:] - colnormB[1, k:]
        minnorm, i = tmp.min(), tmp.argmin()
        q = i + k

        # Column interchange
        if q > k:
            piv[k], piv[q] = piv[q], piv[k]
            colnormB[:, k], colnormB[:, q] = colnormB[:, q], colnormB[:, k]
            B[:, k], B[:, q] = B[:, q].copy(), B[:, k].copy()
        # Compute and apply the Householder transformation  I-tau*v*v'
        if norm(B[k + 1:, k]) > 0:  # otherwise no Householder transformation is needed
            v = B[k:, k].reshape(m - k, 1).copy()  # v = column-vector
            rho = norm(v)  # scalar
            if v[0] >= 0:
                rho = -rho
            v[0] -= rho  # B(k,k)+sgn(B(k,k))*norm(B(k:n,k))
            tao = -1 / (rho * v[0])
            B[k, k] = rho
            B[k:, k + 1:] -= tao * dot(v, dot(v.T, B[k:, k + 1:]))
            # Update y by the Householder transformation
            y[k:] = y[k:, :] - tao * dot(v, dot(v.T, y[k:]))

        # Update colnormB(2,k+1:n)
        colnormB[1, k + 1:] += B[k, k + 1:] * B[k, k + 1:]

    return triu(B[:n, :n]), piv, y


def reduction(B, y):
    """
     [R,Z,y] = reduction(B,y) computes the LLL-QRZ factorization:
               Q'*B*Z = [R; 0] and computes Q'*y. The orthogonal matrix Q
               is not produced. Its goal is to reduce a general integer
               least squares problem to an upper triangular one.

     Input arguments:
        B - m by n real matrix with full column rank
        y - m-dimensional real vector to be transformed to Q'y

     Output arguments:
        R - n by n LLL-reduced upper triangular matrix
        Z - n by n unimodular matrix, i.e., an integer matrix with |det(Z)|=1
        y - m-vector transformed from the input y by Q', i.e., y := Q'*y
    """

    m, n = B.shape
    if m < n:
        raise ValueError("Input matrix is column-rank deficient!")

    m2, n2 = y.shape
    if m != m2 and n2 != 1:
        raise ValueError("Input arguments have a dimension error!")

    # QR with minimum-column pivoting
    R, piv, y_q = qrmcp(B.copy(), y.copy())
    R_q, piv_q, y_q_ = octave.qrmcp(B.copy(), y.copy())
    pprint_two_matrices(R,R_q,name1="R (QR,python)", name2="R (QR, octave)")
    pprint_two_matrices(piv,piv_q,name1="piv (QR,python)", name2="piv (QR, octave)")

    # Obtain the permutation matrix Z
    Z = zeros((n, n))
    for j in xrange(n):
        Z[piv[j]][j] = 1
    # Perform partial LLL reduction on R
    k = 1
    while k < n:
        k1 = k - 1
        zeta = round(R[k1, k] / R[k1, k1])  # NB: division from __future__ is critical here
        alpha = R[k1, k] - zeta * R[k1, k1]

        if R[k1, k1] ** 2 > (1 + 1.e-10) * (alpha ** 2 + R[k, k] ** 2):
            if zeta != 0:
                # Perform a size reduction on R(k-1,k)
                R[k1, k] = alpha
                R[:k-1, k] -= zeta * R[:k-1, k-1]
                Z[:, k] -= zeta * Z[:, k-1]
                # Perform size reductions on R(:k-2,k)
                for i in xrange(k - 2, -1, -1):
                    zeta = round(R[i, k] / R[i, i])
                    if zeta != 0:
                        R[:i, k] = R[:i, k] - zeta * R[:i, i]
                        Z[:, k] = Z[:, k] - zeta * Z[:, i]

            # Permute columns k-1 and k of R
            R[:k+1, k1], R[:k+1, k] = R[:k+1, k].copy(), R[:k+1, k1].copy()
            Z[:, k1], Z[:, k] = Z[:, k].copy(), Z[:, k1].copy()

            # Bring R back to an upper triangular matrix by a Givens rotation
            # TODO: Matlab style, sorry
            r = sqrt(R[k1, k1] ** 2 + R[k, k1] ** 2)
            c = R[k1, k1] / r
            s = R[k, k1] / r
            G = array([[c, s], [-s, c]])
            R[k1, k1] = r
            R[k, k1] = 0
            # Y.K.: tricky code, but note that `k1=k-1`, original code was:
            # R([k1,k],k:n) = G*R([k1,k],k:n);
            R[k1:k + 1, k:] = dot(G, R[k1:k + 1, k:])

            # Apply the Givens rotation to y
            # Y.K.: same here:
            y[k1:k + 1] = dot(G, y[k1:k + 1])

            if k > 1:
                k -= 1
        else:
            k += 1

    return R, Z, y


def search(R, y, p=1):
    """
    z_hat = search(R,y,p) produces p optimal solutions to the upper triangular
           integer least squares problem min_{z}||y-Rz|| by a search algorithm.
    Input arguments:
    R - n by n real nonsingular upper triangular matrix
    y - n-dimensional real vector
    p - the number of optimal solutions with a default value of 1
    Output arguments:
    z_hat - n by p integer matrix (in double precision). Its j-th column
           s the j-th optimal solution, i.e., its residual norm is the j-th
           smallest, so ||y-R*z_hat(:,1)|| <= ...<= ||y-R*z_hat(:,p)||
    """
    m, n = R.shape
    n2, n3 = y.shape
    if m != n or n != n2 or n3 != 1:
        raise ValueError("Input arguments have a matrix dimension error!")

    # Initialization
    z = zeros((n, 1), dtype=int)  # the current point
    c = zeros((n, 1))  # c(k)=(y(k)-R(k,k+1:n)z(k+1:n))/R(k,k)
    d = zeros((n, 1), dtype=int)  # d(k) is used to compute z(k)
    prsd = zeros((n, 1))  # partial squared residual norm for z
    # --> prsd(k)=norm(y(k+1:n)-R(k+1:n,k+1:n)z(k+1:n))^2
    S = zeros((n, n + 1))   # S(:,k) = R(:,k:n)*z(k:n), k=1:n
    z_hat = zeros((n, p), dtype=int)   # the p candidate solutions (or points)
    rsd = zeros((p, 1))     # squared residual norms of the p candidate solutions

    beta = inf  # the initial ellipsoid bound
    ncand = 0   # the initial number of candidate solutions

    nn = n-1
    c[nn] = y[nn] / R[nn, nn]
    z[nn] = int(round(c[nn]))
    gamma = R[nn, nn] * (c[nn] - z[nn])
    if c[nn] > z[nn]:
        d[nn] = 1
    else:
        d[nn] = -1

    k = nn

    while True:
        newprsd = prsd[k] + gamma * gamma
        if newprsd < beta:
            if k != 0:  # move to level k-1
                S[:k+1, k] += R[:k+1, k] * z[k] + S[:k+1, k+1]
                k -= 1
                prsd[k] = newprsd
                c[k] = (y[k] - S[k, k + 1]) / R[k, k]
                z[k] = int(round(c[k]))
                gamma = R[k, k] * (c[k] - z[k])
                if c[k] > z[k]:
                    d[k] = 1
                else:
                    d[k] = -1
            else:  # a new point is found, update the set of candidate solutions
                if ncand < p:  # add the new point
                    z_hat[:, ncand] = z.T
                    rsd[ncand] = newprsd
                    ncand += 1
                    if ncand == p:
                        beta = rsd[-1]
                else:  # insert the new point and remove the worst one
                    i = 0
                    while i < p-1 and rsd[i] < newprsd:
                        i += 1
                    z_hat[:, i:] = [z, z_hat[:, i:-1]]
                    rsd[i:] = [newprsd, rsd[i:-1]]
                    beta = rsd[p-1]
                z[0] = z[0] + d[0]
                gamma = R[0, 0] * (c[0] - z[0])
                if d[0] > 0:
                    d[0] = -d[0] - 1
                else:
                    d[0] = -d[0] + 1
        else:
            if k == n-1:  # the p optimal solutions have been found
                break
            else:  # move back to level k + 1
                k += 1
                z[k] += d[k]
                gamma = R[k, k] * (c[k] - z[k])
                if d[k] > 0:
                    d[k] = -d[k] - 1
                else:
                    d[k] = -d[k] + 1
    return z_hat


def ils1(B, y, p=1):
    m, n = B.shape

    if rank(B) < n:
        raise ValueError("Matrix is rank deficient!")

    # Reduction
    print("=== Reduction ===")
    R, Z_r, y_r = reduction(B.copy(), y.copy())
    _R, _Z_r, _y_r = octave.reduction(B.copy(), y.copy())
    pprint_two_matrices(R, _R,name1="R (python)",name2="R (matlab)")
    pprint_two_matrices(Z_r, _Z_r,name1="Z (python)",name2="Z (matlab)")

    # Search
    _z_hat = search(R.copy(), y_r[:n].copy(), p)
    print("=== Search ===")
    _z_hat_ = octave.search(R.copy(), y_r[:n].copy(), p)

    # Perform the unimodual transformation to obtain the optimal solutions
    return dot(Z_r, _z_hat)


def mils(A, B, y, p=1):
    # x_hat,z_hat = mils(A,B,y,p) produces p pairs of optimal solutions to
    #               the mixed integer least squares problem min_{x,z}||y-Ax-Bz||, 
    #               where x and z are real and integer vectors, respectively.
    #
    # Input arguments:
    #    A - m by k real matrix
    #    B - m by n real matrix
    #          [A,B] has full column rank
    #    y - m-dimensional real vector
    #    p - the number of optimal solutions
    #
    # Output arguments:
    #    x_hat - k by p real matrix
    #    z_hat - n by p integer matrix (in double precision). 
    #           The pair {x_hat(:,j),z_hat(:,j)} is the j-th optimal solution
    #           i.e., its residual is the j-th smallest, so
    #           ||y-A*x_hat(:,1)-B*z_hat(:,1)||<=...<=||y-A*x_hat(:,p)-B*z_hat(:,p)||

    m, k = A.shape
    m2, n = B.shape
    if m != m2 or m != len(y[0]) or len(y[1]) != 1:
        raise ValueError("Input arguments have a matrix dimension error!")

    if rank(A) + rank(B) < k + n:
        raise ValueError("hmmm...")

    Q, R = qr(A)
    Q_A = Q[:, :k]
    Q_Abar = Q[:, k:m]
    R_A = R[:k, :]

    # Compute the p optimal integer least squares solutions
    z_hat = ils(dot(Q_Abar.T, B), dot(Q_Abar.T, y), p)

    # Compute the corresponding real least squares solutions
    # x_hat = solve(R_A, (Q_A.T * (y * ones(p) - B * z_hat)))
    x_hat = lstsq(R_A, dot(Q_A.T, (dot(y, ones(p)) - dot(B, z_hat))))

    return x_hat, z_hat


def pprint_two_matrices(M1, M2, dtype='float', name1='', name2=''):
    """
    Pretty print two matrices (of the same size) in parallel
    """
    if len(M1) != len(M2):
        if len(M1.reshape(1, max(M1.shape))) != len(M2.reshape(1, max(M2.shape))):
            print "Matrices are of different length: %d and %d" % (len(M1), len(M2))
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
    print
    if name1 or name2:
        print "{:^{}}".format(name1, width) + eq + "{:^{}}".format(name2, width)
    for j in xrange(len_M):
        q1, q2 = map(str, _M1[j]), map(str, _M2[j])
        str1 = " ".join(map(lambda x: '{:>{}.{}}'.format(x, max_elem_len + 2, max_elem_len), q1))
        str2 = " ".join(map(lambda x: '{:>{}.{}}'.format(x, max_elem_len + 2, max_elem_len), q2))
        print str1 + ' | ' + str2

if __name__ == "__main__":
    from numpy.random import rand
    from numpy import dot

    m, n = 5, 3
    print "Initital size of matrix B: %d × %d" % (m, n)
    B = rand(m, n)
    # B = array([[-1.38483, 0.53704, 0.14925],  # ans:
    #            [1.05734, 0.61432, 0.94116],  # v1    v2
    #            [-0.33438, -0.13293, -0.60755],  # 1     1
    #            [0.26814, 0.41059, -0.52649],  # -2    -1
    #            [-0.66335, -1.42715, -0.97412]])  # 3     2
    B = array([[ 0.66456276,  0.26352592,  0.853691  ],
               [ 0.91379849,  0.19624231,  0.59176225],
               [ 0.40818463,  0.56173452,  0.84286355],
               [ 0.95138348,  0.5210463,   0.61106613],
               [ 0.68555716,  0.17760152,  0.4324968 ]])
    print "B =\n", B
    # print
    z_true = array([1, -2, 3]).reshape(3, 1)  # column vector
    # y = dot(B, z_true) + 1e-3 * rand(m, 1)
    y = array([-2.01129, 2.65187, -1.89005, -2.13294, -0.73144]).reshape(m, 1)
    # print "y' =", y.T

    p = 2
    Z = ils1(B.copy(), y.copy(), p)

    from oct2py import octave
    from os.path import expanduser as eu
    octave.addpath(eu('~')+'/code/pylgrim/ils')
    Z_oct = octave.ils(B, y, p)

    # print "B =\n", B
    # print "y' =", y.T

    ok = "OK" if (Z == Z_oct).all() else  "WRONG!"

    print "\n=== Solutions: %s ===" % ok

    pprint_two_matrices(Z, Z_oct, name2="From MATLAB:")
