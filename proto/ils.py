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

from numpy import arange, append, array, delete, dot, inf, ones, sign, sqrt, triu, zeros
from numpy.linalg import lstsq, norm, qr
from numpy.linalg import matrix_rank as rank


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
        i = tmp.argmin()
        q = i + k

        # Column interchange
        if q > k:
            piv[k], piv[q] = piv[q].copy(), piv[k].copy()
            colnormB[:, k], colnormB[:, q] = colnormB[:, q].copy(), colnormB[:, k].copy()
            B[:, k], B[:, q] = B[:, q].copy(), B[:, k].copy()
        # Compute and apply the Householder transformation  I-tau*v*v'
        if norm(B[k + 1:, k]) > 0:  # otherwise no Householder transformation is needed
            v = (B[k:, k].copy()).reshape(m - k, 1)  # v = column-vector
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

    # Obtain the permutation matrix Z
    Z = zeros((n, n), dtype=int)
    for j in xrange(n):
        Z[piv[j]][j] = 1
    # Perform partial LLL reduction on R
    k = 1
    while k < n:
        k1 = k - 1
        zeta = int(round(R[k1, k] / R[k1, k1]))  # NB: division from __future__ is critical here
        alpha = R[k1, k] - zeta * R[k1, k1]

        if R[k1, k1] ** 2 > (1 + 1.e-10) * (alpha ** 2 + R[k, k] ** 2):
            if zeta != 0:
                # Perform a size reduction on R(k-1,k)
                R[k1, k] = alpha
                R[:k-1, k] -= zeta * R[:k-1, k-1]
                Z[:, k] -= zeta * Z[:, k-1]
                # Perform size reductions on R(:k-2,k)
                for i in xrange(k - 2, -1, -1):
                    zeta = int(round(R[i, k] / R[i, i]))
                    if zeta != 0:
                        R[:i+1, k] = R[:i+1, k] - zeta * R[:i+1, i]
                        Z[:, k] = Z[:, k] - zeta * Z[:, i]

            # Permute columns k-1 and k of R
            R[:k+1, k1], R[:k+1, k] = R[:k+1, k].copy(), R[:k+1, k1].copy()
            Z[:, k1], Z[:, k] = Z[:, k].copy(), Z[:, k1].copy()

            # Bring R back to an upper triangular matrix by a Givens rotation
            r = sqrt(R[k1, k1] ** 2 + R[k, k1] ** 2)
            c = R[k1, k1] / r
            s = R[k, k1] / r
            G = array([[c, s], [-s, c]])
            R[k1, k1] = r
            R[k, k1] = 0
            R[k1:k + 1, k:] = dot(G, R[k1:k + 1, k:])

            # Apply the Givens rotation to y
            y_q[k1:k + 1] = dot(G, y_q[k1:k + 1])

            if k > 1:
                k -= 1
        else:
            k += 1

    return R, Z, y_q


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
    d[nn] = sign(c[nn] - z[nn])

    k = nn

    while True:
        newprsd = prsd[k] + gamma * gamma
        if newprsd < beta:
            if k != 0:  # move to level k-1
                S[:k+1, k] = R[:k+1, k] * z[k] + S[:k+1, k+1]
                k -= 1
                prsd[k] = newprsd
                c[k] = (y[k] - S[k, k + 1]) / R[k, k]
                z[k] = int(round(c[k]))
                gamma = R[k, k] * (c[k] - z[k])
                d[k] = sign(c[k] - z[k])
            else:  # a new point is found, update the set of candidate solutions
                if ncand < p:  # add the new point
                    z_hat[:, ncand] = z.T
                    rsd[ncand] = newprsd
                    ncand += 1
                    if ncand == p:
                        beta = rsd[-1]
                else:  # insert the new point and remove the worst one
                    # i = 0
                    i = rsd.argmax()
                    # while i < p-1 and rsd[i] <= newprsd:
                    #     i += 1
                    # z_hat[:, i:] = append(z, z_hat[:, :-1], axis=1)
                    z_hat = append(z.copy(), delete(z_hat.copy(), i, axis=1)).reshape((n, p))
                    rsd = append(newprsd, delete(rsd, i))
                    beta = rsd[-1]
                z[0] += d[0]
                gamma = R[0, 0] * (c[0] - z[0])
                d[0] = -d[0] - sign(d[0])
        else:
            if k == n-1:  # the p optimal solutions have been found
                break
            else:  # move back to level k + 1
                k += 1
                z[k] += d[k]
                gamma = R[k, k] * (c[k] - z[k])
                d[k] = -d[k] - sign(d[k])
    return z_hat


def ils(B, y, p=1):
    m, n = B.shape

    if rank(B) < n:
        raise ValueError("Matrix is rank deficient!")

    # Reduction
    R, Z_r, y_r = reduction(B.copy(), y.copy())

    # Search
    _z_hat = search(R.copy(), y_r[:n].copy(), p)
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
    if m != m2 or m != len(y) or len(y[1]) != 1:
        raise ValueError("Input arguments have a matrix dimension error!")

    if rank(A) + rank(B) < k + n:
        raise ValueError("hmmm...")

    Q, R = qr(A, mode='complete')
    Q_A = Q[:, :k]
    Q_Abar = Q[:, k:]
    R_A = R[:k, :]

    # Compute the p optimal integer least squares solutions
    z_hat = ils(dot(Q_Abar.T, B), dot(Q_Abar.T, y), p)

    # Compute the corresponding real least squares solutions
    x_hat = lstsq(R_A, dot(Q_A.T, (dot(y, ones((1, p))) - dot(B, z_hat))))

    return x_hat, z_hat


if __name__ == "__main__":
    from numpy.random import rand

    m, n = 5, 3
    print "Initital size of matrix B: %d × %d" % (m, n)
    B = rand(m, n)
    print "B =\n", B

    z_true = array([1, -2, 3]).reshape(3, 1)  # column vector
    y = dot(B, z_true) + 1e-3 * rand(m, 1)

    print "y' =", y.T

    p = 2
    Z = ils(B.copy(), y.copy(), p)

    ok = "OK" if (Z.T[0] == z_true.T).all() else "WRONG!"
    print "\n=== Solutions: %s ===" % ok
    print Z


    print "MILS:"
    m = 7
    k = 2
    n = 3
    p = 3
    A = rand(m, k)
    B = rand(m, n)
    y = rand(m, 1)

    print "Three pairs of optimal least squares solutions"
    X, Z = mils(A, B, y, p)

    print X
    print Z