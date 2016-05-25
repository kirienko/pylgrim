#!/usr/bin/python
# ! encoding: UTF8

import datetime as dt
import numpy as np
from datetime import datetime
from math import sqrt, sin, cos, atan2
from numpy.numarray import zeros
from gtime import GTime

# Origin of GPS time:
# ori = GTime(year=1980, month=1, day=6,
ori = dt.datetime(year=1980, month=1, day=6,
                  hour=0, minute=0, second=0)

__author__ = 'kirienko'

"""
Prototype of RINEX navigation file parser
"""

mu = 3.986005E14  # WGS84 value of Earth’s universal gravitational parameter [m/s]
c = 2.99792458E8  # GPS value for speed of light [m³/s²]
Omega_dot_e = 7.2921151467 / 1E5  # WGS84 value of Earth’s rotation rate [rad/s]


class Nav(object):
    def __init__(self, data, header):
        """
        Base navigation class that contains raw data, date and PRN number.
        It also has self.A tuple with three elements, but it is important to remember
                that these elements have different meaning for GPS and GLO.

        :param data: lines of observation data: 8 for GPS, 4 for GLO
        :param header: tuple (LEAP, A0, A1) - almanac parameters for time correction
        :return:
        """
        self.raw_data = []
        self.raw_data += data[0][:22].split()  # <-- date
        self.raw_data += [data[0][22 + 19 * i:22 + 19 * (i + 1)] for i in range(3)]
        for d in data[1:]:
            self.raw_data += [d[3 + i * 19:3 + (i + 1) * 19] for i in range(4)]
        self.raw_data = [d.replace('D', 'E') for d in self.raw_data]
        self.A = tuple(map(float, self.raw_data[7:10]))
        self.PRN_number = int(self.raw_data[0])
        self.leap, self.A0, self.A1 = header

        # Time of the observation
        if int(self.raw_data[1]) < 2000: self.raw_data[1] = '20' + self.raw_data[1]
        sec_msec = float(self.raw_data[6])
        self.date = GTime(*(map(int, self.raw_data[1:6]) + [sec_msec]))  # t_oc
        if self.leap is None:
            self.leap = def_leap(self.date)

        try:
            self.eph = map(float, self.raw_data[10:27])  # broadcast ephemeris
            self.week = int(float(self.raw_data[28]))  # GPS week
            epoch = ori + dt.timedelta(days=self.week * 7) #it's a date, i.e. hours = mins = sec = 0
            self.epoch = GTime(*(epoch.timetuple()[:5]))
            seconds = float(self.eph[8])
            self.t_oe = self.epoch + seconds
        except IndexError:
            pass


class NavGPS(Nav):
    def utc2gps(self, t):
        """
        GPS-like time (not to be confused with *GPS week + seconds*)
        :param t:
        :return:
        """
        delta_t = t - self.date
        if abs(delta_t) > 302400:
            if delta_t > 302400:
                while delta_t > 302400:
                    delta_t -= 604800
            else:
                while delta_t < -302400:
                    delta_t += 604800
        # return t + dt.timedelta(seconds= self.A0 + self.A[0] + self.leap + self.A[1]*delta_t)
        # return t + dt.timedelta(seconds=self.A0 + self.leap + self.A1 * delta_t)
        return t - (self.A0 + self.leap + self.A1 * delta_t)
        # return t + dt.timedelta(seconds=self.A0 + self.A1 * delta_t)

    def _time_rel_correction(self, t_sv):
        """
        Relativistic correction of satellite's time (see [1], p.122)
        :param t:
        :return: the relativistic correction term [in seconds]
        """
        # t = self.utc2gps(t_utc)
        t_k = (t_sv - self.epoch) - self.eph[8]  # Time from ephemeris epoch
        return -4.442807622e-10 * self.eph[5] * self.eph[7] * sin(self._ecc_anomaly(t_k))

    def _ecc_anomaly(self, t_k):
        """
        :return: Eccentric anomaly Eₖ(tₖ)
        """
        delta_n = self.eph[2]  # Δn - Mean motion correction [rad/s]
        e = self.eph[5]  # Eccentricity
        sqrt_a = self.eph[7]  # Square root of semimajor axis
        n = sqrt(mu / sqrt_a ** 6) + delta_n  # print " 2) n = sqrt(μ/a³) + Δn = %f [rad/s]" % n
        M_0 = self.eph[3]  # M₀ - Mean anomaly (at time t_oe )
        M_k = M_0 + n * t_k  # print " 4) Mₖ  = M₀  + (n)(tₖ) = %f " % M_k
        E_k = M_k
        for j in xrange(5):
            E_k -= (E_k - e * sin(E_k) - M_k) / (1 - e * cos(E_k))
        # print " 5) Mₖ  = Eₖ  + e sin(Eₖ) = %f " % E_k
        return E_k

    def time_offset(self, t_utc):
        """

        :param t_utc:
        :return: Satellite clock bias
        """
        # t_sv = self.utc2gps(t_utc)
        t_sv = t_utc      # TODO: time means as GPS (not UTC)
        # return self.A[0] + self.A[1]*(t_sv - self.date).total_seconds() + \
        # delta = (t_sv - self.date).total_seconds()
        delta = t_sv - self.date
        # print "δt = %d:%02.f min" % (int(delta/60),abs(delta%int(delta/60)))
        return self.A[0] + self.A[1] * delta + \
               self.A[2] * delta ** 2 + \
               self._time_rel_correction(t_sv)

    def eph2pos(self, t_utc):
        """
        Computes satellite position (ECEF) and clock bias from broadcast ephemeris.
        See [1], Table 4.1 (p. 117) and Table 4.2 (p. 118)
            [2], Appendix E4 (p. 142)
        :param:  t_utc = time of measurement (to be converted to seconds from t_oe)
        :return: r = [Xₖ, Yₖ, Zₖ] - coordinates of satellite in ECEF
        """
        # IODE= self.eph[0]     # Amplitude of sine correction to orbital radius
        C_rs = self.eph[1]  # Amplitude of sine correction to orbital radius
        # delta_n = self.eph[2]     # Δn - Mean motion correction [rad/s]
        # M_0 = self.eph[3]     # M₀ - Mean anomaly (at time t_oe )
        C_uc = self.eph[4]  # Amplitude of cosine correction to argument of latitude
        e = self.eph[5]  # Eccentricity
        C_us = self.eph[6]  # Amplitude of sine correction to argument of latitude
        sqrt_a = self.eph[7]  # Square root of semimajor axis
        t_oe = self.eph[8]  # Reference time of ephemeris from the beginning of GPS week
        C_ic = self.eph[9]  # Amplitude of cosine correction to inclination angle
        Omega = self.eph[10]  # Ω₀ - Longitude of the ascending node (at weekly epoch)
        C_is = self.eph[11]  # Amplitude of sine correction to inclination angle
        i_0 = self.eph[12]  # i₀ - Inclination angle (at time t_oe)
        C_rc = self.eph[13]  # Amplitude of cosine correction to orbital radius
        omega = self.eph[14]  # ω - Argument of perigee (at time t_oe)
        Omega_dot = self.eph[15]  # dΩ/dt - Rate of change of longitude of the ascending node
        IDOT = self.eph[16]  # Rate of change of inclination angle (i.e., di/dt)

        # t = self.utc2gps(t_utc)
        t = t_utc  # TODO: time means as GPS (not UTC)
        a = sqrt_a ** 2  # print " 1) a = (⎷a)² = %.1f [m]" % a
        t_k = t - self.epoch - t_oe  # Time from ephemeris epoch
        # print " 3) tₖ = t - t_oe = %f [s]" % t_k
        E = self._ecc_anomaly(t_k)
        sat_lat_phi = atan2(sqrt(1 - e ** 2) * sin(E), cos(E) - e) + omega
        # print " 7) φ = arctan(α) + ω = %f [rad]" % sat_lat_phi
        d_lat = C_us * sin(2 * sat_lat_phi) + C_uc * cos(2 * sat_lat_phi)
        # print " 8) δφₖ = C_us sin(2φₖ) + C_uc cos(2φₖ) = %f [rad]" % d_lat
        d_r = C_rs * sin(2 * sat_lat_phi) + C_rc * cos(2 * sat_lat_phi)
        # print " 9) δrₖ = C_rs sin(2φₖ) + C_rc cos(2φₖ) = %f [m]" % d_r
        d_i = C_is * sin(2 * sat_lat_phi) + C_ic * cos(2 * sat_lat_phi)
        # print "10) δiₖ = C_is sin(2φₖ) + C_ic cos(2φₖ) = %f [rad]" % d_i
        u_k = sat_lat_phi + d_lat
        # print "11) uₖ = φₖ + δφₖ = %f [rad]" % u_k
        r_k = a * (1 - e * cos(E)) + d_r
        # print "12) rₖ = a (1 − e cos(Eₖ)) + δrₖ = %.1f [m]" % r_k
        i_k = i_0 + IDOT * t_k + d_i
        # print "13) iₖ = i₀ + (di/dt)tₖ  + δiₖ = %f " % i_k
        Omega_k = Omega + (Omega_dot - Omega_dot_e) * t_k - Omega_dot_e * t_oe
        # print "14) Ωₖ = Ω₀ + (Ω̇ - Ω̇ₑ)tₖ  + Ω̇ₑ t₀ₑ = %f " % Omega_k
        x_prime_k = r_k * cos(u_k)  # print "15) x̕ₖ= rₖ cos(uₖ) = %f" % x_prime_k
        y_prime_k = r_k * sin(u_k)  # print "16) y̕ₖ= rₖ sin(uₖ) = %f" % y_prime_k
        x_k = x_prime_k * cos(Omega_k) - y_prime_k * cos(i_k) * sin(Omega_k)
        # print "17) xₖ = x̕ₖcos(Ωₖ) - y̕ₖcos(iₖ)sin(Ωₖ) = %.1f [km]" % (x_k/1000)
        y_k = x_prime_k * sin(Omega_k) + y_prime_k * cos(i_k) * cos(Omega_k)
        # print "18) yₖ = x̕ₖsin(Ωₖ) + y̕ₖcos(iₖ)cos(Ωₖ) = %.1f [km]" % (y_k/1000)
        z_k = y_prime_k * sin(i_k)
        # print "19) zₖ = y̕ₖsin(iₖ) = %.1f [km]" % (z_k/1000)
        r = np.array([x_k, y_k, z_k])  # ECEF coordinates of the satellite in m
        return r


class NavGLO(Nav):
    def __init__(self, data, header):
        Nav.__init__(self, data, header)
        self.TauN, self.GammaN, self.tk = map(float, self.raw_data[7:10])
        y_0, m_0, d_0 = self.date.std.timetuple()[:3]
        self.t_0 = GTime(year=y_0, month=m_0, day=d_0, hour=0, minute=0)
        self.r = np.array([self.eph[0], self.eph[4], self.eph[8]]) * 1e3  # R₀ = (X₀, Y₀, Z₀)
        self.v = np.array([self.eph[1], self.eph[5], self.eph[9]]) * 1e3  # V₀ = (V_x₀, V_y₀, V_z₀)
        self.acc = np.array([self.eph[2], self.eph[6], self.eph[10]]) * 1e3  # a₀ = (a_x₀, a_y₀, a_z₀)
        self.t_oe = self.t_0 + self.tk

        self.RE_GLO = 6378136.0  # radius of earth [m]               ref [2]
        self.MU_GLO = 3.986004418e14  # gravitational constant [m³/s²]         ref [2]
        self.J2_GLO = 1.08262575e-3  # 2nd zonal harmonic of geopot   ref [2]
        self.OMGE_GLO = 7.292115e-5  # earth angular velocity [rad/s] ref [2]
        self.ERREPH_GLO = 5.0  # error of glonass ephemeris [m]

    def deq(self, x):
        """
        Differential equation solution
        :param x: array of length 6 that contains coordinates and velocities
        :return:
        """

        r2 = np.dot(x[:3], x[:3])
        r3 = r2 * sqrt(r2)
        omg2 = self.OMGE_GLO * self.OMGE_GLO

        if r2 <= 0:
            return zeros(6)

        deq_a = 1.5 * self.J2_GLO * self.MU_GLO * self.RE_GLO**2 / r2 / r3  # /* 3/2*J2*mu*Ae^2/r^5 */
        deq_b = 5.0 * x[2] * x[2] / r2  # /* 5*z^2/r^2 */
        deq_c = -self.MU_GLO / r3 - deq_a * (1.0 - deq_b)  # /* -mu/r^3-a(1-b) */
        xdot0_2 = x[3:6]
        xdot_3 = (deq_c + omg2) * x[0] + 2.0 * self.OMGE_GLO * x[4] + self.acc[0]
        xdot_4 = (deq_c + omg2) * x[1] - 2.0 * self.OMGE_GLO * x[3] + self.acc[1]
        xdot_5 = (deq_c - 2.0 * deq_a) * x[2] + self.acc[2]
        return np.append(xdot0_2, np.array([xdot_3, xdot_4, xdot_5]))

    def eph2pos(self, time):
        """
        ECEF coordinates of the satellite (based algorithm described in RTKlib manual)
        :param time: GTime instance
        """
        TSTEP = 60.0  # integration step glonass ephemeris (s)

        t = (time-self.t_0) - self.tk

        if t > 1800:
            print "Warning: probably wrong interpolation period"

        # glo_xyzt = np.append(r_0, -TauN + GammaN * t)  # coordinates in ECEF and clock bias

        x = np.append(self.r, self.v)   # x = (\vec{r},\vec{v})
        tt = np.sign(t) * TSTEP
        i = 0
        # print "time = ", t, "\tt_0 =", self.t_0, "\tt_k =", self.tk, "\tdate =", self.date, "\ttt =", tt
        while abs(t) > 1E-9:
            if abs(t) < TSTEP:
                tt = t
            k1 = self.deq(x)
            w = x + k1 * tt / 2
            k2 = self.deq(w)
            w = x + k2 * tt / 2
            k3 = self.deq(w)
            w = x + k3 * tt
            k4 = self.deq(w)
            x += (k1 + 2*k2 + 2*k3 + k4) * tt/6
            t -= tt
            i += 1
        # return np.append(x[:3], -self.TauN + self.GammaN * t)  # coordinates in ECEF and clock bias
        return x[:3]        # coordinates in ECEF and clock bias

    def time_offset(self, time):
        t = (time-self.t_0).total_seconds() - self.tk
        return -(-self.TauN + self.GammaN * t)


class PreciseNav:
    def __init__(self, date, sat_position):
        self.date = date
        self.t_oe = self.date
        self.xyzt = np.array(map(float, sat_position))  # [km, km, km, mcs]

    def eph2pos(self, time=0):
        """
        :param time: unnecessary parameter
        :return: ECEF coordinates of the satelite
        """
        return self.xyzt[:3] * 1e3

    def time_offset(self, time=0):
        """
        :param time: unnecessary parameter
        :return: satellite clock bias [sec]
        """
        return self.xyzt[3] / 1e6


def def_leap(date):
    """
    Return the number of leap seconds since 6/i/1980
    :param date: datetime instance
    :return: leap seconds for date (int)
    """
    if date < datetime(1981, 6, 30, 23, 59, 59):
        return 0
    leap_list = [(1981, 6, 30), (1982, 6, 30), (1983, 6, 30),
                 (1985, 6, 30), (1987, 12, 31), (1989, 12, 31),
                 (1990, 12, 31), (1992, 6, 30), (1993, 6, 30),
                 (1994, 6, 30), (1995, 12, 31), (1997, 6, 30),
                 (1998, 12, 31), (2005, 12, 31), (2008, 12, 31),
                 (2012, 6, 30), (2015, 6, 30)]
    leap_dates = map(lambda x: datetime(x[0], x[1], x[2], 23, 59, 59), leap_list)
    for j in xrange(len(leap_dates[:-1])):
        if leap_dates[j] < date < leap_dates[j + 1]:
            return j + 1
    return len(leap_dates)
