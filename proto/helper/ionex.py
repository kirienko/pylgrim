#! encoding: utf8
from proto.parse_rinex import get_header_body, get_header_line
import datetime as dt
import numpy as np


def get_int_from_header(hdr, seq):
    """
    Returns the first int from the line that contains `seq` of lines `hdr`.
    In fact, _header_ here may not be header of RINEX/IONEX, just some set of lines.
    """
    return int(get_header_line(hdr, seq).split()[0])


def closest_in_list(lst, val, num = 2):
    """
    Returns two (`num` in general) closest values of `val` in list `lst`
    """
    return sorted(lst, key=lambda x: abs(x - val))[:num]

class IonexMap:
    def __init__(self, exp, data):
        self.exp = 10 ** (exp)
        self.grid_TEC = np.array([], dtype='uint16')
        date = map(int, data[0].split()[:6])
        self.date = dt.datetime(*date)
        self.lats = np.array([])
        for j, line in enumerate(data[1:]):
            if "LAT" in line:
                lat, lon1, lon2, dlon, h = map(float,[line[x:x+6] for x in xrange(2,32,6)])
                self.lats = np.append(self.lats, lat)
                row_length = (lon2 - lon1)/dlon + 1
                # next_lines_with_numbers = int(np.ceil(row_length / 16))
                next_lines_with_numbers = int(row_length / 16)
                last_line_len = int(row_length % 16)
                row = np.array([], dtype='int16')
                # print j, lat, lon1, lon2, dlon, h
                for i in xrange(next_lines_with_numbers):
                    row = np.append(row, np.array(map(int, [data[j+2+i][5*x:5*x+5] for x in xrange(16)]),dtype='int16'))
                if last_line_len:
                    row = np.append(row, np.array(map(int, [data[j+2 + next_lines_with_numbers][5*x:5*x+5] for x in xrange(last_line_len)]), dtype='int16'))
                if len(self.grid_TEC) > 0:
                    self.grid_TEC = np.vstack((self.grid_TEC, row))
                else:
                    self.grid_TEC = np.append(self.grid_TEC, row)
        self.lons = np.linspace(lon1, lon2, row_length)

    @staticmethod
    def find_nearest(lst, val):
        return (np.abs(lst - val)).argmin()

    def get_TEC(self, pos):
        """
        Returns TEC in a position `pos` (lat, lon) of ionosphere
        :param pos:
        :return:
        """
        lat = self.find_nearest(self.lats, pos[0])
        lon = self.find_nearest(self.lats, pos[1])
        return self.grid_TEC[lat][lon] * self.exp

    @staticmethod
    def round_to_grid(number, base):
        return int(base * round(float(number)/base))


def parse_ionex(ionex_file):
    """
    :param ionex_file: path to the IONEX file
    :return: list of `IonexMap` objects,
    """
    header, body = get_header_body(ionex_file)

    exponent = get_int_from_header(header, "EXPONENT")
    maps_count = get_int_from_header(header, "MAPS IN FILE")

    # =============
    # Separate maps
    # =============
    map_start_idx = []
    map_end_idx = []

    for j, line in enumerate(body):
        if "START OF TEC MAP" in line:
            map_start_idx += [j]
        elif "END OF TEC MAP" in line:
            map_end_idx += [j]
    if maps_count != len(map_start_idx):
        raise LookupError("Parsing error: the number of maps in the header "
                          "is not equal to the number of maps in the body.")
    if len(map_start_idx) != len(map_end_idx):
        raise IndexError("Starts end ends numbers are not equal.")
    # print map_start_idx, len(map_start_idx)
    map_dates = []
    for i in xrange(maps_count):
        date = map(int, body[map_start_idx[i] + 1].split()[:6])
        map_dates += [dt.datetime(*date)]
    print map_dates
    # def get_TEC(time):
    #     if time in map_dates:
    #         return map_dates.index(time)
    #     elif time > map_dates[0] and time < map_dates[-1]:
    #         # find two maps closest in time
    #         closest = closest_in_list(map_dates, time)
    #         return map(lambda x: map_dates.index(x),closest)
    #     else:
    #         raise IndexError("No IONEX map for this time: %s" % str(time))
    maps = []
    for m in xrange(maps_count):
        iono_map = body[map_start_idx[m]+1:map_end_idx[m]]
        maps += [IonexMap(exponent, iono_map)]
    return maps


if __name__ == "__main__":
    from os.path import expanduser

    ionex = expanduser("~") + "/code/pylgrim/test_data/iter_0/igsg0010.16i"
    M = parse_ionex(ionex)
