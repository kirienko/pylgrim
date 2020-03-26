from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import map
from builtins import range
import re
from collections import defaultdict
from math import ceil
from .gtime import GTime
from .nav_data import NavGPS, NavGLO, PreciseNav
from .obs_data import ObsGPS
from proto.helper.parsing_utils import get_header_line, get_header_body

__author__ = 'kirienko'


def parse_rinex(path):
    """
    Parse RINEX-file and returns a timeline, i.e. array of tuples
        (time_of_event, Object_instance),
    where Object_instance is either Nav or Observ depending on the type of
    RINEX file (detected automatically).
    :param path: path to RINEX file
    :return: timeline
    """

    header, body = get_header_body(path)

    # Define RINEX file properties
    rinex_version = re.findall('[0-9]' + '\.?' + '[0-9]*', header[0][:60])[0]
    print("\n\nParsing %s:" % path)
    print("RINEX version:", rinex_version)
    pat_file = re.compile("nav|obs", re.IGNORECASE)
    pat_satel = re.compile("gps|glo|mix", re.IGNORECASE)
    rinex_type = re.findall(pat_file, header[0][:60])[0].lower()
    print("RINEX type:", rinex_type)
    try:
        satel_type = re.findall(pat_satel, header[0][:60])[0].lower()
    except IndexError:
        pass
    if float(rinex_version) < 2.11:
        print("Warning: Rinex version is too old.")
        satel_type = 'gps'

    if rinex_type == 'nav':
        # Define UTC conversion data
        utc = get_header_line(header, 'delta-utc')
        try:
            LEAP = int(get_header_line(header, 'leap')[:60])
            A0, A1 = [float(h.replace('D', 'E')) for h in [utc[:22], utc[22:41]]]
            T, W = list(map(int, utc[42:60].split()))
            print("Satellites:", satel_type)
            print("utc: A0 = %e, A1 = %e, T = %d, W = %d, leap seconds: %d" % (A0, A1, T, W, LEAP))
        except TypeError:
            A0, A1 = 0., 0.
            if satel_type == 'mix':
                LEAP = None
            else:
                LEAP = 0

    prefixes = {'gps': 'G', 'glo': 'R', 'mix': 'M'}
    sat_prefix = prefixes[satel_type]

    if rinex_type == 'nav' and satel_type != 'mix':

        if satel_type == 'gps':
            if len(body) % 8 != 0:
                print("Warning: wrong length of NAV file")
            nav_list = [NavGPS(body[i * 8:(i + 1) * 8], (LEAP, A0, A1)) for i in range(len(body) // 8)]
        elif satel_type == 'glo':
            if len(body) % 4 != 0:
                print("Warning: wrong length of NAV file")
            nav_list = [NavGLO(body[i * 4:(i + 1) * 4], (LEAP, A0, A1)) for i in range(len(body) // 4)]
        nav_dict = defaultdict(list)
        for obj in nav_list:
            nav_dict[sat_prefix + "%02d" % obj.PRN_number] += [obj]
        return nav_dict

    elif rinex_type == 'obs':
        obs_types = get_header_line(header, "TYPES OF OBSERV").split()
        number_of_obs_types = int(obs_types[0])
        lpo = int(ceil(float(number_of_obs_types) / 5))  # lines per one observation
        obs_types = obs_types[1:1 + number_of_obs_types]

        observations = []
        for j, h in enumerate(body):
            if ('G' in h or 'R' in h) and h[31] in '0123456789':
                satellite_count = int(h[30:32])
                if satellite_count > 12:  # sometimes it happens!
                    observations += [ObsGPS(body[j:j + satellite_count * lpo + 2], obs_types)]
                else:
                    observations += [ObsGPS(body[j:j + satellite_count * lpo + 1], obs_types)]
        return observations
    else:
        raise NotImplementedError


def parse_sp3(path):
    """
    Parse SP3-file that contains precise ephemeris
    :param path: path to *.sp3-file
    :return: array of Nav-objects
    """
    print("\nParsing %s:" % path)
    with open(path) as fd:
        data = fd.readlines()
    nav_dict = defaultdict(list)
    for j, d in enumerate(data):
        if d[0] == '*':
            split = d.split()[1:]
            y, m, d, H, M = list(map(int, split[:-1]))
            s = float(split[-1])
            date = GTime(y, m, d, H, M, s)
        elif d[0] == 'P' and date:  # GPS satellites
            prn, x, y, z, t = d[2:].split()[:5]
            nav_dict[d[1] + "%02d" % int(prn)] += [PreciseNav(date, (x, y, z, t))]
        else:
            continue
    return nav_dict
