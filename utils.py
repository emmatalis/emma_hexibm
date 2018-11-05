import os
import random
import bisect as _bisect
import itertools
import math
import numpy as np
import sys

def flatten(l):
    return [item for sublist in l for item in sublist]


class Random36(random.Random):
    "Show the code included in the Python 3.6 version of the Random class"

    def choices(self, population, weights=None, *, cum_weights=None, k=1):
        """Return a k sized list of population elements chosen with replacement.
        If the relative weights or cumulative weights are not specified,
        the selections are made with equal probability.
        """
        random = self.random
        if cum_weights is None:
            if weights is None:
                _int = int
                total = len(population)
                return [population[_int(random() * total)] for i in range(k)]
            cum_weights = list(itertools.accumulate(weights))
        elif weights is not None:
            raise TypeError('Cannot specify both weights and cumulative weights')
        if len(cum_weights) != len(population):
            raise ValueError('The number of weights does not match the population')
        bisect = _bisect.bisect
        total = cum_weights[-1]
        return [population[bisect(cum_weights, random() * total)] for i in range(k)]


random36 = Random36()


def select_weighted(d, lam=lambda x: x):
    """
    weighted random sample from dictionary
    :param d: dict
    :param lam: callable - modify dict value to weighting
    :return: obj, dictionary value
    """
    prob_sum = int(sum([lam(x) for x in d.values()]) - 1)
    offset = random.randint(0, prob_sum)
    for key, value in d.items():
        if offset < lam(value):
            return key
        offset -= lam(value)


def convert_size(size_bytes):
    if size_bytes == 0:
        return "0B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return "%s %s" % (s, size_name[i])


class HexNeighbors:
    """ Provides methods to get the indicies of neighboring cells in hexagonal grid
    Methods for finding neighbors in a hexgaonal grid. Assumes hexgonal grid uses
    axial coordinate system x,y,z (http://www.redblobgames.com/grids/hexagons/#coordinates).
    Coordinates are calculated for a rings of hexagons, with caching to speed the process
    for repeat calls to the class.
    """

    def __init__(self):
        self.cache = {}

    def get_ring_offset(self, distance):
        """
        Return the indices for a ring of hexagons at 'distance' from an origin hexagon of (0,0,0)
        :param distance: int
        :return: list
        """
        if distance in self.cache:
            return self.cache[distance]
        else:
            coords_positive = list(zip(range(0, distance + 1), range(distance, -1, -1)))
            coords_negative = list(zip(range(-distance, 1), range(0, -distance - 1, -1)))

            all_coords = list(set(itertools.chain([(x, y, -distance) for (x, y) in coords_positive],
                                                  [(-distance, y, z) for (z, y) in coords_positive],
                                                  [(x, -distance, z) for (z, x) in coords_positive],
                                                  [(x, y, distance) for (x, y) in coords_negative],
                                                  [(x, distance, z) for (x, z) in coords_negative],
                                                  [(distance, y, z) for (z, y) in coords_negative])))
            self.cache[distance] = all_coords
            print("cache extended to {} units ({})".format(distance, convert_size(sys.getsizeof(self.cache))))

            return all_coords

    def get_radius_offset(self, end_distance, start_distance=0):
        """
        Return indices of all hexagons within radius 'distance'.
        :param distance: int
        :return: list
        """
        return flatten([self.get_ring_offset(i) for i in range(start_distance, end_distance + 1)])

    def get_ring_coords(self, distance, origin=(0, 0, 0)):
        """
        Return indices of all hexagons in a ring at 'distance' from specified origin.
        :param distance: int
        :param origin: tuple
        :return: list
        """
        x_, y_, z_ = origin
        return [(x_ + x, y_ + y, z_ + z) for x, y, z in self.get_ring_offset(distance)]

    def get_radius_coords(self, end_distance, origin=(0, 0, 0), start_distance=0):
        x_, y_, z_ = origin
        return [(x_ + x, y_ + y, z_ + z) for x, y, z in
                self.get_radius_offset(end_distance=end_distance, start_distance=start_distance)]


class IncrementOutputFile:
    """
    Returns a numbered file name that increments with each call
    Note: By default does not check if file exists and is not thread
    safe
    """

    def __init__(self, base: str, ndigits: int = 3, ext: str = '', sep: str = '_', warn: bool = False):
        """
        :param base:
        :param ndigits:
        :param ext:
        :param sep:
        :param warn:
        """
        self.base = base
        self.ndigits = ndigits
        self.count = 0
        self.warn = warn
        self.sep = sep
        self.ext = ext

    def __call__(self):
        """
        Creates uniquely numbered file, increments count, and returns path to new
        file.
        :return: string
        """
        path = "{0}{1}{2:0{3}d}{4}".format(self.base, self.sep, self.count, self.ndigits, self.ext)
        if self.warn and os.path.exists(path):
            raise IOError("File already exists")
        self.count += 1
        return path


class UniqueOutputDirectory:
    def __init__(self, root):
        dirs = os.listdir(root)
        self.root = root
        self.max_count = 0
        for directory in dirs:
            try:
                self.max_count = max([int(directory), self.max_count])
            except:
                self.max_count = 0

    def __call__(self, *args, **kwargs):
        self.max_count += 1
        path = os.path.join(self.root, str(self.max_count))
        if not os.path.isdir(path):
            os.mkdir(path)
        return path


def inv_logit(p):
    """Return inverse logit of p"""
    return np.exp(p) / (1 + np.exp(p))


def notify(message, title=" "):
    """
    Uses ubuntu notify to issue system notification
    :param message: string
    :param title: string
    :return: None
    """
    os.system("notify-send {} {}".format(title, message))
    return None

def estBetaParams(mu, var):
    """
    Moment matching for Beta distribution
    :param mu:
    :param var:
    :return:
    """
    alpha = ((1 - mu) / var - 1 / mu) * mu ** 2
    beta = alpha * (1 / mu - 1)
    return alpha, beta