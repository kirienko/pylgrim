#!/usr/bin/python
#! encoding: UTF8
__author__ = 'kirienko'

from nav_data import Nav
from obs_data import Observ
import re
from collections import defaultdict

def parse_rinex(path):
    """
    Parse RINEX-file and returns a timeline, i.e. array of tuples
        (time_of_event, Object_instance),
    where Object_instance is either Nav or Observ depending on the type of RINEX file (detected automatically).
    :param path: path to RINEX file
    :return: timeline
    """

    with open(path) as fd:
        data = fd.readlines()
        for j,d in enumerate(data):
            if "END OF HEADER" in d:
                header_end = j
                break
    header, body = data[:header_end], data[header_end+1:]

    # Define RINEX file properties
    rinex_version = re.findall('[0-9]'+'\.'+'[0-9]+',header[0][:60])[0]
    pat_file = re.compile("nav|obs",re.IGNORECASE)
    pat_satel = re.compile("gps|glo|mix",re.IGNORECASE)
    rinex_type = re.findall(pat_file,header[0][:60])[0].lower()
    satel_type = re.findall(pat_satel,header[0][:60])[0].lower()
    print "RINEX version:", rinex_version
    print "RINEX type:", rinex_type
    print "Satellites:", satel_type

    prefixes = {'gps':'G', 'glo':'R'}
    sat_prefix = prefixes[satel_type]

    if rinex_type == 'nav' and satel_type != 'mix': #FIXME: it should work for all satellite types
        if len(body) % 8 != 0:
            print "Warning: wrong length of NAV file"
        nav_list = [Nav(body[i*8:(i+1)*8]) for i in xrange(len(body)/8)]
        nav_dict = defaultdict(list)
        for obj in nav_list:
            nav_dict[sat_prefix+"%02d"%obj.PRN_number] += [obj]
        return nav_dict
    else:
        raise NotImplementedError


if __name__ == "__main__":

    navigations = parse_rinex('../test_data/test.n')
    for k,v in navigations.items(): print k, [str(vv.date)for vv in v]

    g16 = navigations['G16']
    z2,z1 = g16[1],g16[2]
    t1,t2 = z1.eph[7], z2.eph[7]
    print z1.date, z2.date