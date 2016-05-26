#! encoding: utf8
import datetime as dt
import numpy as np
from parsing_utils import get_header_body, get_int_from_header

nu_1_G, nu_2_G = 1575.42e6, 1227.6e6


def closest_in_list(lst, val, num=2):
    """
    Returns two (`num` in general) closest values of `val` in list `lst`
    """
    idxs = sorted(lst, key=lambda x: abs(x - val))[:num]
    return sorted(map(lambda x: list(lst).index(x), idxs))


class IonexMap:
    def __init__(self, exp, data):
        self.exp = 40.30E16 / (nu_1_G ** 2) * 10 ** (exp)
        self.grid_TEC = np.array([], dtype='uint16')
        date = map(int, data[0].split()[:6])
        self.date = dt.datetime(*date)
        self.lats = np.array([])
        for j, line in enumerate(data[1:]):
            if "LAT" in line:
                lat, lon1, lon2, dlon, h = map(float, [line[x:x + 6] for x in xrange(2, 32, 6)])
                self.lats = np.append(self.lats, lat)
                row_length = (lon2 - lon1) / dlon + 1  # total number of values of longitudes
                next_lines_with_numbers = int(np.ceil(row_length / 16))
                elems_in_row = [min(16, int(row_length - i * 16)) for i in range(next_lines_with_numbers)]
                row = np.array([], dtype='int16')
                for i, elem in enumerate(elems_in_row):
                    row = np.append(row, np.array(map(int, [data[j + 2 + i][5 * x:5 * x + 5] for x in xrange(elem)]),
                                                  dtype='int16'))
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
        Returns TEC in a position `pos`  of ionosphere
        :param pos: (lat, lon) [deg, deg]
        :return:
        """
        if pos[0] in self.lats and pos[1] in self.lons:
            lat = self.find_nearest(self.lats, pos[0])
            lon = self.find_nearest(self.lons, pos[1])
            E = self.grid_TEC[lat][lon] * self.exp
            return E
        lat_idxs = closest_in_list(self.lats, pos[0])
        lon_idxs = closest_in_list(self.lons, pos[1])
        lat0, lat1 = self.lats[lat_idxs[0]], self.lats[lat_idxs[1]]
        lon0, lon1 = self.lons[lon_idxs[0]], self.lons[lon_idxs[1]]
        dlat = lat1 - lat0
        dlon = lon1 - lon0
        p = float(pos[0] - lat0) / dlat
        q = float(pos[1] - lon0) / dlon
        (E00, E10), (E01, E11) = self.grid_TEC[lat_idxs[0]:lat_idxs[1] + 1, lon_idxs[0]:lon_idxs[1] + 1]
        return ((1 - p) * (1 - q) * E00 + p * (1 - q) * E01 + (1 - p) * q * E10 + p * q * E11) * self.exp

    @staticmethod
    def round_to_grid(number, base):
        return int(base * round(float(number) / base))


def parse_ionex(ionex_file):
    """
    :param ionex_file: path to the IONEX file
    :return: TEC interpolation function `f( (lat,lon), datetime )`
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
    map_dates = []
    for i in xrange(maps_count):
        date = map(int, body[map_start_idx[i] + 1].split()[:6])
        map_dates += [dt.datetime(*date)]

    maps = []
    for m in xrange(maps_count):
        iono_map = body[map_start_idx[m] + 1:map_end_idx[m]]
        maps += [IonexMap(exponent, iono_map)]

    def interpolate_maps(position, time):
        """
        Interpolate TEC at position at the moment `time`
        :param position: (lat,lon) [deg]
        :param time: datetime instance
        :return:
        """
        if time in map_dates:
            return maps[map_dates.index(time)].get_TEC(position)
        elif time not in map_dates and time > map_dates[0] and time < map_dates[-1]:
            # find two maps closest in time
            closest = closest_in_list(map_dates, time)
            m0, m1 = maps[closest[0]], maps[closest[1]]
            t0, t1 = m0.date, m1.date
            dt0, dt1 = (time - t0).total_seconds(), (time - t1).total_seconds()
            theta = dt0 / (t1 - t0).total_seconds()
            e0, e1 = m0.get_TEC([position[0], position[1] + dt0 / 240]), m1.get_TEC(
                [position[0], position[1] + dt1 / 240])
            return (1 - theta) * e0 + theta * e1
        else:
            raise IndexError("No IONEX map for this time: %s" % str(time))

    return interpolate_maps
