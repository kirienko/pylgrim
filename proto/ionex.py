#! encoding: utf8
from numpy import array, append
def round_to_grid(number, base):
    return int(base * round(float(number)/base))


def find_VMF_coeffs(ah_file, aw_file, coords):
    if ah_file.replace('ah', '') != aw_file.replace('aw', ''):
        print "WARNING: file names do not match:", ah_file, aw_file
    crds = coords[:]    # just a copy
    coeffs = []
    for f in (ah_file, aw_file):
        with open(f) as fd:
            data = fd.readlines()
            values = array([])
            for line in data[1:]:
                values = append(values, array([float(x) for x in line.strip().split()]))
        # print len(values)
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
        grid_coords = [round_to_grid(x, grid[j]) for j, x in enumerate(crds)]
        coeffs += [values[int((lat_range[0] - grid_coords[0]) / grid[0])][int(grid_coords[1] / grid[1])]/10E6]

    return coeffs

if __name__ == "__main__":
    ahf = '../test_data/iter_0/VMF_ah16001.h00'
    awf = '../test_data/iter_0/VMF_aw16001.h00'
    # ahf = '../test_data/iter_0/ah08079.h00'   # ah = 0.00121328
    # awf = '../test_data/iter_0/aw08079.h00'   # aw = 0.00043331
    print find_VMF_coeffs(ahf, awf, (46,15))
