#! encoding: utf8

import datetime as dt


class GTime(object):
    """
    Time class based on datetime, but with seconds of type ``float``
    """

    def __init__(self, year=None,
                 month=None,
                 day=None,
                 hour=None,
                 minute=None,
                 second=0.0):
        self.std = dt.datetime(year=year, month=month, day=day, hour=hour,
                               minute=minute)
        self.sec = float(second)

    def __sub__(self, other):
        """
        Equivalent of dt.timedelta.total_seconds()
        """
        return (self.std - other.std).total_seconds() + self.sec - other.sec

    def __add__(self, other):
        """
        If we add a float, treat it like seconds
        """
        if isinstance(other, float) or isinstance(other, int):
            totalsec = self.sec + other
            s_std = self.std + dt.timedelta(minutes=int(totalsec / 60))
            s_sec = totalsec % 60
            s_std = list(s_std.timetuple()[:5]) + [s_sec]
            summed = GTime(*s_std)
            return summed
        else:
            print("Type of addend: " + str(type(other)))
            raise NotImplementedError

    def __lt__(self, other):
        if self.std < other.std:
            return True
        elif self.std > other.std:
            return False
        else:
            return self.sec < other.sec

    def __str__(self):  # I'm ugly and I know it
        return str(self.std)[:-2] + "%02d" % int(self.sec) + ("%.9f" % (self.sec - int(self.sec)))[1:]
