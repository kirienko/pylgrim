#! encoding: utf8
from proto.parse_rinex import get_header_body, get_header_line
import datetime as dt


def get_int_from_header(hdr, seq):
    """
    Returns the first int from the line that contains `seq` of lines `hdr`.
    In fact, _header_ here may not be header of RINEX/IONEX, just some set of lines.
    """
    return int(get_header_line(hdr, seq).split()[0])


if __name__ == "__main__":
    from os.path import expanduser

    ionex_file = expanduser("~") + "/code/pylgrim/test_data/iter_0/igsg0010.16i"

    header, body = get_header_body(ionex_file)

    exponent = get_int_from_header(header, "EXPONENT")
    interval = get_int_from_header(header, "INTERVAL")
    maps_count = get_int_from_header(header, "MAPS IN FILE")
    print interval, exponent

    # =============
    # Separate maps
    # =============
    map_start_idx = []

    for j, line in enumerate(body):
        if "START OF TEC MAP" in line:
            map_start_idx += [j]
    if maps_count != len(map_start_idx):
        raise LookupError("Parsing error: the number of maps in the header "
                          "is not equal to the number of maps in the body.")
    print map_start_idx, len(map_start_idx)
    map_dates = []
    for i in xrange(maps_count):
        date = map(int, body[map_start_idx[i] + 1].split()[:6])
        map_dates += [dt.datetime(*date)]
    print map_dates
