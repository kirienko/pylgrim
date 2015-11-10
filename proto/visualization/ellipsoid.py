#! encoding: UTF8
import numpy as np
from numpy import sin,cos, tan, pi, ones, zeros

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from coord.ecef import ecef_to_lat_lon_alt


def satellites(pos, sat_pos, sat_names=''):
    a = 6378137.0
    b = 6356752.314245

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_aspect('equal',anchor='C')
    # ax.set_axis_off()

    u = np.linspace(0, 2 * pi, 100)
    v = np.linspace(0, pi, 100)

    x = a * np.outer(cos(u), sin(v))
    y = a * np.outer(sin(u), sin(v))
    z = b * np.outer(ones(np.size(u)), cos(v))
    ax.plot_surface(x, y, z, rstride=5, cstride=2, color='blue',
                    alpha=0.5,
                    linewidth=0.05)
    ax.plot(a * cos(u),
            a * sin(u), 0, '-r',linewidth = 0.5)  # <-- equator
    ax.plot(a * cos(v-pi/2), zeros(100), b * sin(v-pi/2), '-r',linewidth=0.5)

    xx,yy,zz = [[p[i] for p in sat_pos] for i in range(3)]
    ax.scatter3D(xx,yy,zz,color='k',s=10)

    if sat_names:
        for i in range(len(sat_names)):
            ax.text(xx[i],yy[i],zz[i],sat_names[i])
    if not isinstance(pos,list):
        pos = list(pos)
    ax.scatter3D(*pos,color='r',s=10)

    # cone
    elev_mask = np.deg2rad(60)
    perpendicular = [(0,5*pos[i]) for i in range(3)]
    ax.plot(*perpendicular)

    # rotate cone on phi and theta

    def rotate(X,Y,alpha):
        R = np.matrix([[cos(alpha),-sin(alpha)],
                       [sin(alpha),cos(alpha)]])
        print len(X),len(Y)
        M = np.matrix([X, Y])
        X1, Y1 = R * M
        # print X1,Y1
        return X1.getA()[0], Y1.getA()[0]

    phi, theta, h = ecef_to_lat_lon_alt(pos,deg=False)

    cone_x = np.linspace(0,4*a,100)
    cone_y = cone_x
    cone_z = cone_x
    X,Y1 = rotate(cone_x,cone_y,theta)
    # X1,Z1 = rotate(X,cone_z,phi)

    # cone_y = cone_x * np.outer(cos(u), 1/tan(elev_mask))
    cone_z = cone_x * np.outer(sin(u), 1/tan(elev_mask))
    cone_y = X * np.outer(cos(u), 1/tan(elev_mask))
    cone_z = X * np.outer(sin(u), 1/tan(elev_mask))

    ax.plot_surface(cone_x, cone_y, cone_z,
                    rstride=5, cstride=5, color='y',
                    alpha=0.2, linewidth=0.01)
    ax.set_xlim3d(-5*a,5*a)
    ax.set_ylim3d(-5*a,5*a)
    ax.set_zlim3d(-5*a,5*a)
    plt.show()