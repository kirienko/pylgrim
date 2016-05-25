#! encoding: utf8
from numpy import array, append, cos, sin, deg2rad, pi
from datetime import datetime
"""
 - - - - - - - - - 
  V M F 1 _ H T 
 - - - - - - - - - 

 This routine is part of the International Earth Rotation and 
 Reference Systems Service (IERS) Conventions software collection. 

 This subroutine determines the Vienna Mapping Function 1 (VMF1) (Boehm et al. 2006). 

    :------------------------------------------: 
    :                                          : 
    :                 IMPORTANT                : 
    :                                          : 
    :  This version uses height correction!    : 
    :  It has to be used with the VMF Grid     : 
    :  located at the website mentioned in     : 
    :  the Notes.                              : 
    :__________________________________________: 

 In general, Class 1, 2, and 3 models represent physical effects that 
 act on geodetic parameters while canonical models provide lower-level 
 representations or basic computations that are used by Class 1, 2, or 
 3 models. 

 Status: Class 1 model 

    Class 1 models are those recommended to be used a priori in the 
    reduction of raw space geodetic data in order to determine 
    geodetic parameter estimates. 
    Class 2 models are those that eliminate an observational 
    singularity and are purely conventional in nature. 
    Class 3 models are those that are not required as either Class 
    1 or 2. 
    Canonical models are accepted as is and cannot be classified as a 
    Class 1, 2, or 3 model. 

 Given: 
    AH             d      Hydrostatic coefficient a (Note 1) 
    AW             d      Wet coefficient a (Note 1) 
    DMJD           d      Modified Julian Date 
    DLAT           d      Latitude given in radians (North Latitude) 
    HT             d      Ellipsoidal height given in meters 
    ZD             d      Zenith distance in radians 

 Returned: 
    VMF1H          d      Hydrostatic mapping function (Note 2) 
    VMF1W          d      Wet mapping function (Note 2) 

 Notes: 

 1) The coefficients can be obtained from the primary website 
    http://ggosatm.hg.tuwien.ac.at/DELAY/ or the back-up website 
    http://www.hg.tuwien.ac.at/~ecmwf1/. 

 2) The mapping functions are dimensionless scale factors. 

 Test case: 
    given input: AH   = 0.00127683D0 
                 AW   = 0.00060955D0 
                 DMJD = 55055D0 
                 DLAT = 0.6708665767D0 radians (NRAO, Green Bank, WV) 
                 HT   = 824.17D0 meters 
                 ZD   = 1.278564131D0 radians 

    expected output: VMF1H = 3.423513691014495652D0 
                     VMF1W = 3.449100942061193553D0 

 References: 

    Boehm, J., Werl, B., and Schuh, H., (2006), 
    "Troposhere mapping functions for GPS and very long baseline 
    interferometry from European Centre for Medium-Range Weather 
    Forecasts operational analysis data," J. Geophy. Res., Vol. 111, 
    B02406, doi:10.1029/2005JB003629 

    Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), 
    IERS Technical Note No. 36, BKG (2010) 

 Revisions: 
 2005 October 02 J. Boehm     Original code 
 2009 August 17 B.E. Stetzler Added header and copyright 
 2009 August 17 B.E. Stetzler More modifications and defined twopi 
 2009 August 17 B.E. Stetzler Provided test case 
 2009 August 17 B.E. Stetzler Capitalized all variables for FORTRAN 77 
                              compatibility 
 2010 September 08 B.E. Stetzler   Provided new primary website to obtain 
                                   VMF coefficients 
--------------------------------------------------------------------- 
--------------------------------------------------------------------- 
    Reference day is 28 January 1980 
    This is taken from Niell (1996) to be consistent 
--------------------------------------------------------------------- 
"""

def round_to_grid(number, base):
    return int(base * round(float(number)/base))


def find_VMF_coeffs(ah_file, aw_file, coords):
    """

    :param ah_file:
    :param aw_file:
    :param coords: (lat, lon) in degrees; it's OK if coords = (lat,lon,hgt)
    :return:
    """
    if ah_file.replace("ah", '') != aw_file.replace('aw', ''):
        print "WARNING: file names do not match:", ah_file, aw_file
    crds = coords[:]    # just a copy
    coeffs = []
    for f in (ah_file, aw_file):
        with open(f) as fd:
            data = fd.readlines()
            values = array([])
            for line in data[1:]:
                values = append(values, array([float(x) for x in line.strip().split()]))
        header = map(float, data[0].split())
        lat_range = header[0], header[1]
        lon_range = header[2], header[3]
        grid = header[4], header[5]
        lat_len = abs(lat_range[0] - lat_range[1])/grid[0] + 1
        lon_len = abs(lon_range[0] - lon_range[1])/grid[1] + 1
        values = values.reshape((lat_len, lon_len))

        if abs(lat_range[0] - lat_range[1]) < 180 or \
            abs(lon_range[0] - lon_range[1]) < 360:
                print "WARNING: ionospheric map does not cover all the Earth!"
        while crds[1] < 0:
            crds[1] += 360
        grid_coords = [round_to_grid(x, grid[j]) for j, x in enumerate(crds[:2])]
        coeffs += [values[int((lat_range[0] - grid_coords[0]) / grid[0])][int(grid_coords[1] / grid[1])]/10E7]

    return tuple(coeffs)


def cont_fraction(el, a, b, c):
    """

    :param el: elevation angle [rad]
    :param a, b, c: coefficients
    :return: continued fraction
    """
    sine = sin(el)
    beta = b / (sine + c)
    gamma = a / (sine + beta)
    topcon = a / (b / (c + 1.) + 1.) + 1.
    return topcon / (sine + gamma)


def vmf(pos, time, elev, coeffs):
    """
    Vienna mapping function (aka VMF1_HT)
    :param pos: receiver position in (lat, lon, alt) [rad, rad, m]
    :param time: time of observation, datetime object
    :param elev: elevation angle [rad]
    :param coeffs: hydrostatic (ah) and wet (aw) coefficients, e.g.: (0.00121328, 0.00043331)
    :return: tropospheric delay in meters
    """
    ah, aw = coeffs
    ep = datetime(2000, 1, 1, 12, 0, 0)
    mjd = 51544.5 + (time - ep).total_seconds() / 86400.0
    doy = mjd - 44239. - 28
    bh = .0029
    c0h = .062
    if pos[0] < 0.:
        # southern hemisphere
        psi = pi
        c11h = .007
        c10h = .002
    else:
        # northern hemisphere
        psi = 0.
        c11h = .005
        c10h = .001

    ch = c0h + ((cos(doy / 365.25 * 6.283185307179586476925287 + psi) + 1.) *
                c11h / 2. + c10h) * (1. - cos(pos[0]))
    vmf1h = cont_fraction(elev, ah, bh, ch)
    # Compute the height correction
    #   (Niell, A. E. "Global mapping functions for the atmosphere delay at radio wavelengths."
    #   Journal of Geophysical Research: Solid Earth 101.B2 (1996): 3227-3246.)
    #   DOI: 10.1029/95JB03048
    a_ht__ = 2.53e-5
    b_ht__ = .00549
    c_ht__ = .00114
    hs_km__ = float(pos[2]) / 1e3
    ht_corr_coef__ = 1. / sin(elev) - cont_fraction(elev, a_ht__, b_ht__, c_ht__)
    ht_corr__ = ht_corr_coef__ * hs_km__
    vmf1h += ht_corr__
    bw = 0.00146    # taken from Böhm, Schuh, Atmospheric effects in space geodesy, 2013, page 103
    cw = 0.04391    # taken from Böhm, Schuh, Atmospheric effects in space geodesy, 2013, page 103
    vmf1w = cont_fraction(elev, aw, bw, cw)

    """
    +---------------------------------------------------------------------- 

     Copyright (C) 2008 
     IERS Conventions Center 

     ================================== 
     IERS Conventions Software License 
     ================================== 

     NOTICE TO USER: 

     BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS 
     WHICH APPLY TO ITS USE. 

     1. The Software is provided by the IERS Conventions Center ("the 
        Center"). 

     2. Permission is granted to anyone to use the Software for any 
        purpose, including commercial applications, free of charge, 
        subject to the conditions and restrictions listed below. 

     3. You (the user) may adapt the Software and its algorithms for your 
        own purposes and you may distribute the resulting "derived work" 
        to others, provided that the derived work complies with the 
        following requirements: 

        a) Your work shall be clearly identified so that it cannot be 
           mistaken for IERS Conventions software and that it has been 
           neither distributed by nor endorsed by the Center. 

        b) Your work (including source code) must contain descriptions of 
           how the derived work is based upon and/or differs from the 
           original Software. 

        c) The name(s) of all modified routine(s) that you distribute 
           shall be changed. 

        d) The origin of the IERS Conventions components of your derived 
           work must not be misrepresented; you must not claim that you 
           wrote the original Software. 

        e) The source code must be included for all routine(s) that you 
           distribute.  This notice must be reproduced intact in any 
           source distribution. 

     4. In any published work produced by the user and which includes 
        results achieved by using the Software, you shall acknowledge 
        that the Software was used in obtaining those results. 

     5. The Software is provided to the user "as is" and the Center makes 
        no warranty as to its use or performance.   The Center does not 
        and cannot warrant the performance or results which the user may 
        obtain by using the Software.  The Center makes no warranties, 
        express or implied, as to non-infringement of third party rights, 
        merchantability, or fitness for any particular purpose.  In no 
        event will the Center be liable to the user for any consequential, 
        incidental, or special damages, including any lost profits or lost 
        savings, even if a Center representative has been advised of such 
        damages, or for any claim by any third party. 

     Correspondence concerning IERS Conventions software should be 
     addressed as follows: 

                        Gerard Petit 
        Internet email: gpetit[at]bipm.org 
        Postal address: IERS Conventions Center 
                        Time, frequency and gravimetry section, BIPM 
                        Pavillon de Breteuil 
                        92312 Sevres  FRANCE 

        or 

                        Brian Luzum 
        Internet email: brian.luzum[at]usno.navy.mil 
        Postal address: IERS Conventions Center 
                        Earth Orientation Department 
                        3450 Massachusetts Ave, NW 
                        Washington, DC 20392 


    ----------------------------------------------------------------------- 
    """
    return vmf1h, vmf1w
