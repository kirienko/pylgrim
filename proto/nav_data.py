#!/usr/bin/python
#! encoding: UTF8

from time import mktime, time
import datetime as dt
from math import sqrt, sin, cos, atan2

__author__ = 'kirienko'


"""
Prototype of RINEX navigation file parser


References:
[1] Mohinder S. Grewal, Angus P. Andrews, Chris G. Bartone,
    Global Navigation Satellite Systems Inertial Navigation and Integration,
    Wiley-Interscience, 2013
    http://ru.bookzz.org/book/2086274/def272

"""

def get_gps_week_number(now = None):
    """
    :returns GPS week number from UNIX time, or datetime timetuple;
        if argument is None returns the number of current week
    """
    unix_to_gps = 315964800
    if isinstance(now, dt.datetime):
        ## suppose datetime timetuple
        now = int(mktime(now.timetuple()))
    elif isinstance(now,int) and now > unix_to_gps:
        ## suppose timestamp
        pass
        #return (now - unix_to_gps)/(60*60*24*7)
    elif now is None:
        now = int(time())
        #return (now - unix_to_gps)/(60*60*24*7)
    else:
        print "\n\tError: cannot determine GPS week"
        exit(1)
    return (now - unix_to_gps)/(60*60*24*7)

mu = 3.986005E14                # WGS84 value of Earth’s universal gravitational parameter [m/s]
c = 2.99792458E8                # GPS value for speed of light [m³/s²]
Omega_dot_e = 7.2921151467/1E5  # WGS84 value of Earth’s rotation rate [rad/s]

class Nav():
    def __init__(self, data):
        self.raw_data = []
        self.raw_data += data[0][:22].split()   # <-- date
        self.raw_data +=[data[0][22+19*i:22+19*(i+1)] for i in range(3)]
        for d in data[1:]:
            self.raw_data += [d[3+i*19:3+(i+1)*19] for i in range(4)]
        self.raw_data = [d.replace('D','E') for d in self.raw_data]
        self.PRN_number = int(self.raw_data[0])

        # Time of the observation
        if int(self.raw_data[1]) < 2000: self.raw_data[1] = '20' + self.raw_data[1]
        sec_msec = "%.3f" % float(self.raw_data[6])
        s,ms = map(int,sec_msec.split('.'))
        self.date = dt.datetime(*(map(int,self.raw_data[1:6])+[s,ms]))
        # self.week = get_gps_week_number(self.date)    # do we need this?

        self.A   = map(float,self.raw_data[7:10])
        self.eph = map(float,self.raw_data[11:27])  # broadcast ephemeris
        self.epoch = self.date - dt.timedelta(seconds=self.eph[7]) # beginning of current GPS week
        self.now = self.date - self.epoch

    def eph2pos(self):
        """
        Computes satellite position (ECEF) and clock bias from broadcast ephemeris.\
        See [1], Table 4.1 (p. 117) and Table 4.2 (p. 118)
        :return:
        """
        C_rs      = self.eph[0]     # Amplitude of sine correction to orbital radius
        delta_n   = self.eph[1]     # Δn - Mean motion correction [rad/s]
        M_0       = self.eph[2]     # M₀ - Mean anomaly (at time t_0e )
        C_uc      = self.eph[3]     # Amplitude of cosine correction to argument of latitude
        e         = self.eph[4]     # Eccentricity
        C_us      = self.eph[5]     # Amplitude of sine correction to argument of latitude
        sqrt_a    = self.eph[6]     # Square root of semimajor axis
        t_0e      = self.eph[7]     # Reference time of ephemeris from the beginning of GPS week
        C_ic      = self.eph[8]     # Amplitude of cosine correction to inclination angle
        Omega     = self.eph[9]     # Ω₀ - Longitude of the ascending node (at weekly epoch)
        C_is      = self.eph[10]    # Amplitude of sine correction to inclination angle
        i_0       = self.eph[11]    # i₀ - Inclination angle (at time t_0e)
        C_rc      = self.eph[12]    # Amplitude of cosine correction to orbital radius
        omega     = self.eph[13]    # ω - Argument of perigee (at time t_0e)
        Omega_dot = self.eph[14]    # dΩ/dt - Rate of change of longitude of the ascending node
        IDOT      = self.eph[15]    # Rate of change of inclination angle (i.e., di/dt)

        a = sqrt_a**2           # print " 1) a = (⎷a)² = %.1f [m]" % a
        n = sqrt(mu/a**3) + delta_n  # print " 2) n = sqrt(μ/a³) + Δn = %f [rad/s]" % n
        t_k = self.now - t_0e          # Time from ephemeris epoch
        # print " 3) tₖ = t - t_0e = %f [s]" % t_k.seconds
        M_k = M_0 + n * t_k     # print " 4) Mₖ  = M₀  + (n)(tₖ) = %f " % M_k
        E = M_k
        for j in xrange(5):
            E -= (E - e * sin(E) - M_k) / (1 - e * cos(E))
        # print " 5) Mₖ  = Eₖ  + e sin(Eₖ) = %f " % E
        sat_lat_phi = atan2(sqrt(1 - e**2)*sin(E),cos(E)-e) + omega
        # print " 7) φ = arctan(α) + ω = %f [rad]" % sat_lat_phi
        d_lat = C_us*sin(2*sat_lat_phi) + C_uc*cos(2*sat_lat_phi)
        # print " 8) δφₖ = C_us sin(2φₖ) + C_uc cos(2φₖ) = %f [rad]" % d_lat
        d_r = C_rs*sin(2*sat_lat_phi) + C_rc*cos(2*sat_lat_phi)
        # print " 9) δrₖ = C_rs sin(2φₖ) + C_rc cos(2φₖ) = %f [m]" % d_r
        d_i = C_is*sin(2*sat_lat_phi) + C_ic*cos(2*sat_lat_phi)
        # print "10) δiₖ = C_is sin(2φₖ) + C_ic cos(2φₖ) = %f [rad]" % d_i
        u_k = sat_lat_phi + d_lat
        # print "11) uₖ = φₖ + δφₖ = %f [rad]" % u_k
        r_k = a*(1 - e * cos(E)) + d_r
        # print "12) rₖ = a (1 − e cos(Eₖ)) + δrₖ = %.1f [m]" % r_k
        i_k = i_0 + IDOT * t_k + d_i
        # print "13) iₖ = i₀ + (di/dt)tₖ  + δiₖ = %f " % i_k
        Omega_k = Omega + (Omega_dot-Omega_dot_e)*t_k + Omega_dot_e*t_0e #TODO
        # print "14) Ωₖ = Ω₀ + (Ω̇ - Ω̇ₑ)tₖ  + Ω̇ₑ t₀ₑ = %f " % Omega_k
        x_prime_k = r_k * cos(u_k)  # print "15) x̕ₖ= rₖ cos(uₖ) = %f" % x_prime_k
        y_prime_k = r_k * sin(u_k)  # print "16) y̕ₖ= rₖ sin(uₖ) = %f" % y_prime_k
        x_k = x_prime_k * cos(Omega_k) - y_prime_k * cos(i_k) * sin(Omega_k)
        # print "17) xₖ = x̕ₖcos(Ωₖ) - y̕ₖcos(iₖ)sin(Ωₖ) = %.1f [km]" % (x_k/1000)
        y_k = x_prime_k * sin(Omega_k) - y_prime_k * cos(i_k) * cos(Omega_k)
        # print "18) yₖ = x̕ₖsin(Ωₖ) + y̕ₖcos(iₖ)cos(Ωₖ) = %.1f [km]" % (y_k/1000)
        z_k = y_prime_k * sin(i_k)
        # print "19) zₖ = y̕ₖsin(iₖ) = %.1f [km]" % (z_k/1000)
        r = [x_k, y_k, z_k]     # ECEF coordinates of the satellite

        return r

if __name__ == "__main__":
    with open('../test_data/test.n') as fd:
        data = fd.readlines()

    header_end_marker = "END OF HEADER"
    for j,d in enumerate(data):
        if header_end_marker in d:
            header_end = j
    header, body = data[:header_end], data[header_end+1:]

    for i in xrange(len(body[:8])):
        print body[i],
        if (i+1) % 8 == 0: print

    o = Nav(body[:8])
    print o.eph2pos(10.0)