import re


def get_header_line(headr, proprty):
    """
    :param headr: the header of the RINEX-file
    :param proprty: string-like property to search for (e.g. 'delta-utc')
    :return: the string of the ``headr`` containing ``property``
    """
    pattern = re.compile(proprty, re.IGNORECASE)
    for d in headr:
        if pattern.search(d):
            return d


def get_header_body(file_path):
    """
    Opens `file_path`, reads file and returns header and body
    separated with "END OF HEADER"
    :param file_path: path to RINEX-like file
    :return: header, body (arrays of lines)
    """
    with open(file_path) as fd:
        data = fd.readlines()
        for j, d in enumerate(data):
            if "END OF HEADER" in d:
                header_end = j
                break
    return data[:header_end], data[header_end + 1:]


def get_int_from_header(hdr, seq):
    """
    Returns the first int from the line that contains `seq` of lines `hdr`.
    In fact, _header_ here may not be header of RINEX/IONEX, just some set of lines.
    """
    return int(get_header_line(hdr, seq).split()[0])