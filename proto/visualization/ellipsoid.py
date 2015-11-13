#! encoding: UTF8
import numpy as np
from numpy import sin,cos, tan, pi, ones, zeros

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from coord.ecef import ecef_to_lat_lon_alt, sat_elev


def satellites(pos, sat_pos, sat_names=''):
    '''
    Visualize satellites above the rover
    :param pos: rover's position (ecef)
    :param sat_pos: array of satellites positions (each position in ecef)
    :param sat_names: [optional] satellites names
    :return: None
    '''
    a = 6378137.0       # Major semi-axis
    b = 6356752.314245  # Minor semi-axis
    elev_mask = np.deg2rad(30)  # elevation mask
    num = 100           # number of points to render
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_aspect('equal',anchor='C')
    ax.set_axis_off()

    u = np.linspace(0, 2 * pi, num)
    v = np.linspace(0, pi, num)

    x = a * np.outer(cos(u), sin(v))
    y = a * np.outer(sin(u), sin(v))
    z = b * np.outer(ones(num), cos(v))
    ax.plot_surface(x, y, z, rstride=5, cstride=2, color='blue',
                    alpha=0.5,linewidth=0.05)
    ax.plot(a * cos(u),
            a * sin(u), 0, '-r',linewidth = 0.5)  # <-- equator
    ax.plot(a * cos(v-pi/2), zeros(num), b * sin(v-pi/2), '-r',linewidth=0.5)

    xx,yy,zz = [[p[i] for p in sat_pos] for i in range(3)]

    for idx, sat in enumerate(sat_pos):
        if sat_elev(pos,sat,deg=False) > elev_mask:
            ax.scatter3D(*sat,color='k',s=10) # <-- satellites
        else:
            ax.scatter3D(*sat, alpha=0.5,
                         color='gray',s=10) # <-- satellites
            # print "Bad satellite:",sat, sat_elev(pos,sat)

    if sat_names:
        for i in range(len(sat_names)):
            ax.text(xx[i],yy[i],zz[i],sat_names[i])
    if not isinstance(pos,list):
        pos = list(pos)
    ax.scatter3D(*pos[:3],color='r',s=10)

    perpendicular = [(pos[i],3.5*pos[i]) for i in range(3)]
    ax.plot(*perpendicular,color='y')

    def rotate(X,Y,Z,phi,theta):
        '''
        Euler's rotations
        :param X:
        :param Y:
        :return:
        '''
        alpha = theta
        beta  = pi/2
        gamma = phi
        D = np.matrix([[cos(alpha),-sin(alpha), 0],
                       [sin(alpha), cos(alpha), 0],
                       [    0,          0,      1]])
        C= np.matrix([[cos(beta), 0,-sin(beta)],
                       [0,     1,       0],
                       [sin(beta),0, cos(beta)]])
        B = np.matrix([[cos(gamma), sin(gamma), 0],
                       [-sin(gamma),cos(gamma), 0],
                       [    0,          0,      1]])

        X1, Y1, Z1 = [],[],[]
        for i in xrange(num):
            vec = np.array([X[i],Y[i],Z[i]])
            XYZ = B * C * D * vec
            X1.extend(XYZ[0].tolist())
            Y1.extend(XYZ[1].tolist())
            Z1.extend(XYZ[2].tolist())
        return map(np.array,[X1, Y1, Z1])

    phi, theta, h = ecef_to_lat_lon_alt(pos,deg=False)

    # cone:
    t = np.linspace(0,2.5*a,num)      # just a parameter
    cone_x = np.outer(ones(num),t)
    cone_y = t * np.outer(cos(u), 1/tan(elev_mask))
    cone_z = t * np.outer(sin(u), 1/tan(elev_mask))
    # rotate cone on phi and theta
    Xx,Yy, Zz = rotate(cone_x,cone_y,cone_z,phi,theta)

    ax.plot_surface(Xx+pos[0], Yy+pos[1], Zz+pos[2],
                    rstride=5, cstride=5, color='y',
                    alpha=0.1, linewidth=0.01)
    ax.set_xlim3d(-4*a,4*a)
    ax.set_ylim3d(-4*a,4*a)
    ax.set_zlim3d(-4*a,4*a)
    plt.show()