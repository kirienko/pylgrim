#! encoding: UTF8

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.basemap import Basemap

def draw_inset(plt,m,pos,lat,lon):
    ax = plt.subplot(111)
    zoom = 50
    axins = zoomed_inset_axes(ax, zoom, loc=1)
    m.plot(lon,lat,'.b--',zorder=10,latlon=True)
    m.scatter(lon,lat,      # longitude first!
              latlon=True,  # lat and long in degrees
              zorder=11)   # on top of all
    x1, y1 = m(lon[1]-0.005,lat[0]-0.0025)
    x2, y2 = m(lon[1]+0.005,lat[0]+0.0025)
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)

    plt.xticks(visible=False)
    plt.yticks(visible=False)
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

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

    map_width = map_height = 1e6    # optimized for the scale: 1e5

    m = Basemap(width=map_width,height=map_height,\
                resolution='c',     # 'c', 'l', 'i', 'h', 'f'
                projection=proj,lon_0=pos[1],lat_0=pos[0])
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
    # draw parallels and meridians.
    m.drawparallels(np.arange(int(pos[0])-5,int(pos[0])+5,1.),labels=[1,0,0,1])
    m.drawmeridians(np.arange(int(pos[1])-5,int(pos[1])+5,1.),labels=[1,0,0,1])
    m.drawmapboundary(fill_color='aqua')
    m.plot(lon,lat,'.b--',zorder=10,latlon=True)

    m.scatter(lon,lat,      # longitude first!
              color='r',
              latlon=True,  # lat and long in degrees
              zorder=11)   # on top of all
    x0,y0 = m(pos[1],pos[0])    # center os the plot
    lon_lab, lat_lab = m(x0+0.4*map_width,y0-0.4*map_width,inverse=True)
    m.drawmapscale(
        lon_lab,lat_lab,  # where to place scale
        pos[1],pos[0],          # where to measure scale
        map_width/1e4,                     # length
        units='km', fontsize=10,
        barstyle='fancy', labelstyle='simple',
        fillcolor1='w', fillcolor2='#000000',
        fontcolor='#000000',
        zorder=20)
    m.drawmapscale(pos[1],pos[0],pos[0],pos[1],20,barstyle='simple',zorder=50)
    plt.title(u'Map centered at (%.1fN, %.1fE)' % pos)
    # draw inset with starting position
    # draw_inset(plt,m,pos,lat,lon)
    plt.show()