#! encoding: UTF8

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def on_map(positions,proj='cass'):
    '''
    The supported projections are:
     cea              Cylindrical Equal Area
     mbtfpq           McBryde-Thomas Flat-Polar Quartic
     aeqd             Azimuthal Equidistant
     sinu             Sinusoidal
     poly             Polyconic
     omerc            Oblique Mercator
     gnom             Gnomonic
     moll             Mollweide
     lcc              Lambert Conformal
     tmerc            Transverse Mercator
     nplaea           North-Polar Lambert Azimuthal
     gall             Gall Stereographic Cylindrical
     npaeqd           North-Polar Azimuthal Equidistant
     mill             Miller Cylindrical
     merc             Mercator
     stere            Stereographic
     eqdc             Equidistant Conic
     rotpole          Rotated Pole
     cyl              Cylindrical Equidistant
     npstere          North-Polar Stereographic
     spstere          South-Polar Stereographic
     hammer           Hammer
     geos             Geostationary
     nsper            Near-Sided Perspective
     eck4             Eckert IV
     aea              Albers Equal Area
     kav7             Kavrayskiy VII
     spaeqd           South-Polar Azimuthal Equidistant
     ortho            Orthographic
     cass             Cassini-Soldner
     vandg            van der Grinten
     laea             Lambert Azimuthal Equal Area
     splaea           South-Polar Lambert Azimuthal
     robin            Robinson
    :param pos:
    :param proj:
    :return:
    '''
    # llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
    # are the lat/lon values of the lower left and upper right corners
    # of the map.
    # resolution = 'i' means use intermediate resolution coastlines.
    # lon_0, lat_0 are the central longitude and latitude of the projection.
    # m = Basemap(llcrnrlon=-10.5,llcrnrlat=49.5,urcrnrlon=3.5,urcrnrlat=59.5,
    #             resolution='i',projection=proj,lon_0=-4.36,lat_0=54.7)
    # can get the identical map this way (by specifying width and
    # height instead of lat/lon corners)
    lat, lon = [p[0] for p in positions], [p[1] for p in positions]
    pos = (sum(lat)/len(lat), sum(lon)/len(lon))

    map_width = map_height = 1e5

    m = Basemap(width=map_width,height=map_height,\
               resolution='i',projection=proj,lon_0=pos[1],lat_0=pos[0])
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
    # draw parallels and meridians.
    m.drawparallels(np.arange(int(pos[0])-5,int(pos[0])+5,1.),labels=[1,0,0,1])
    m.drawmeridians(np.arange(int(pos[1])-5,int(pos[1])+5,1.),labels=[1,0,0,1])
    m.drawmapboundary(fill_color='aqua')
    # print pos.tolist()
    m.plot(lon,lat,'.b--',zorder=10,latlon=True)

    m.scatter(lon,lat,      # longitude first!
              color='r',
              latlon=True,  # lat and long in degrees
              zorder=11)   # on top of all

    m.drawmapscale(
        pos[1]+0.8,pos[0]-0.4,  # where to place scale
        pos[1],pos[0],          # where to measure scale
        10,                     # length
        units='km', fontsize=10,
        barstyle='fancy', labelstyle='simple',
        fillcolor1='w', fillcolor2='#000000',
        fontcolor='#000000',
        zorder=20)
    m.drawmapscale(pos[1],pos[0],pos[0],pos[1],20,barstyle='simple',zorder=50)
    plt.show()