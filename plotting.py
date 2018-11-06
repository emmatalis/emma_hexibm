import csv
from ast import literal_eval as make_tuple
import fiona
import numpy as np
import itertools
import random
import pickle
import rasterio
from rasterio import features
from shapely.geometry import shape
import subprocess
import re
import math
import sys


def zeropad_filename(fn,n=4):
    return re.sub(
        '(\d+)',
        lambda m: m.group().zfill(n),
        fn)

flatten = lambda l: [item for sublist in l for item in sublist]



def convert_size(size_bytes):
    if size_bytes == 0:
        return "0B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return "%s %s" % (s, size_name[i])


class HexNeighbors:
    """ Provides cached methods to get the indicies of neighboring cells in hexagonal grid

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

    def get_radius_offset(self, distance):
        """
        Return indices of all hexagons within radius 'distance'.
        :param distance: int
        :return: list
        """
        return flatten([self.get_ring_offset(i) for i in range(distance + 1)])

    def get_ring_coords(self, distance, origin=(0, 0, 0)):
        """
        Return indices of all hexagons in a ring at 'distance' from specified origin.
        :param distance: int
        :param origin: tuple
        :return: list
        """
        x_, y_, z_ = origin
        return [(x_ + x, y_ + y, z_ + z) for x, y, z in self.get_ring_offset(distance)]

    def get_radius_coords(self, distance, origin=(0, 0, 0)):
        x_, y_, z_ = origin
        return [(x_ + x, y_ + y, z_ + z) for x, y, z in self.get_radius_offset(distance)]


def inv_logit(p):
    """Return inverse logit of p"""
    return np.exp(p) / (1 + np.exp(p))


class IBMmap:
    def __init__(self):
        self.hexes = {}
        self.neighbor_manager = HexNeighbors()
        # Sim annealing variable
        self.T = 1
        self.T_min = 0.00001
        self.alpha = 0.9

    def pickle(self, path):
        with open(path, 'wb') as hex_pickle:
            return pickle.dump(self.hexes, hex_pickle)

    def from_pickle(self, path):
        with open(path, 'rb') as hex_pickle:
            self.hexes = pickle.load(hex_pickle)

    def get_neighbors(self, hex, distance):
        return [self.hexes.get(hex) for hex in self.neighbor_manager.get_radius_coords(distance, hex.axial_coords)]

    def get_neighbors_ring(self, hex, distance):
        return [self.hexes.get(hex) for hex in self.neighbor_manager.get_ring_coords(distance, hex.axial_coords)]

    def get_suitability(self, hexes=None):
        if not hexes:
            self.suitability = np.nansum(
                [self.hexes[hex].get_quality() for hex in self.hexes if self.hexes[hex].properties['occupied'] == 1])
            return self.suitability
        else:
            return np.nansum(
                [self.hexes[hex].get_quality() for hex in hexes if self.hexes[hex].properties['occupied'] == 1])

    def get_occupied(self):
        self.occupied = set([hex for hex in self.hexes if self.hexes[hex].properties['occupied'] == 1])

    def get_unoccupied(self):
        self.unoccupied = set([hex for hex in self.hexes if self.hexes[hex].properties['occupied'] == 0])

    def set_keys_list(self):
        self.keys_list = list(self.hexes)

    def accept(self, new, old):
        return np.exp((new - old) / self.T) > random.random()

    def test_accept(self):
        backupT = self.T
        self.T = 1
        assert self.accept(1, 0) == True
        print("annealing acceptance test 1 passed")

        self.T = 0.0000001
        assert self.accept(1, 0) == True
        print("annealing acceptance test 2 passed")

        self.T = 10
        assert self.accept(0, 1) == True
        print("annealing acceptance test 3 passed")

        self.T = 0.0000001
        assert self.accept(0, 1) == False
        print("all annealing acceptance tests passed")
        self.T = backupT

    def update_T(self):
        self.T = max(self.T * self.alpha, self.T_min)

    def switch(self):
        found_occupied = False
        found_unoccupied = False
        while not found_occupied:
            occ = random.sample(self.keys_list, 1)[0]
            found_occupied = self.hexes[occ].properties['occupied'] == 1
        while not found_unoccupied:
            unocc = random.sample(self.keys_list, 1)[0]
            found_unoccupied = self.hexes[unocc].properties['occupied'] == 0
        previous_suitability = self.hexes[occ].get_quality()
        new_suitability = self.hexes[unocc].get_quality()
        # if new_suitability >= previous_suitability:
        #     self.hexes[occ].properties['occupied'] = 0
        #     self.hexes[unocc].properties['occupied'] = 1
        if self.accept(new_suitability, previous_suitability):
            self.hexes[occ].properties['occupied'] = 0
            self.hexes[unocc].properties['occupied'] = 1


class Hex:
    def __init__(self, grid, axial, properties):
        self.grid = grid
        self.axial_coords = axial
        self.properties = properties
        self.fon = False

    def __repr__(self):
        return "Hex at {} {} {}".format(*self.axial_coords)

    @property
    def fono(self):
        if not self.fon:
            self.fon = self.grid.get_neighbors(self, 1)
        return np.sum([x.properties['occupied'] for x in self.fon if x is not None])

    def get_neighbors(self, distance):
        return self.grid.get_neighbors(self, distance)

    def get_neighbors_ring(self, distance):
        return self.grid.get_neighbors_ring(self, distance)

    def get_quality(self):
        p = -3.6 - 0.0004 * self.properties['flow'] - 0.0005 * self.properties['distance'] + 0.0005 * self.properties[
            'elevation'] + 1.6 * self.fono
        return inv_logit(p)


def load_state(path):
    with open(path, 'r') as f:
        reader = csv.reader(f)
        hexkeys = [make_tuple(key[0]) for key in reader]
    return hexkeys


def parse_filename(self, filename: str):
    return int(filename.split('_')[0])

class CstLineHolder:
    def __init__(self):
        self.properties = {}
        with fiona.collection('coastline/beagle_cst.shp', 'r') as layer:
            for element in layer:
                self.properties['geom'] =  shape(element['geometry'])
        self.properties['col'] = 100

def load_ages(path):
    with open(path, 'r') as f:
        reader = csv.reader(f)
        ages = [int(row[1]) for row in reader]
    return ages

class Plotter:
    def __init__(self):
        self.transform = [-2398667.8006973956,
                          0.06537704951100207,
                          0.0,
                          1700808.706766952,
                          0.0,
                          -0.3837805797099703]
        self.shape = (3735, 18220)

        self.keys = None
        self.template = IBMmap()
        self.template.from_pickle('hex_map_sim_anneal')
        self.template.hexes['COASTLINE'] = CstLineHolder()
        print("template loaded.")

        self.temp_files = []

    def _shapes(self):
        return ((self.template.hexes[key].properties['geom'], age) for key,age in zip(self.keys,self.ages))

    def rasterise_keys(self, keys):
        self.keys = ['COASTLINE']
        self.keys += keys
        return features.rasterize(self._shapes(), out_shape=self.shape, transform=self.transform)

    def reconstruct(self, path):
        outpath = os.path.split(path)[0]
        print("outpath : ",outpath)
        filename = os.path.splitext(os.path.split(path)[1])[0]
        out_filename = zeropad_filename(filename)
        print("outfilename : ",out_filename)
        out_filename = '{}.tif'.format(os.path.join(outpath,out_filename))
        print("outfilename : ",out_filename)
        keys = load_state(path=path)
        self.ages = [255] + load_ages(path=path)
        image = self.rasterise_keys(keys)
        print('rasterised')
        with rasterio.open(out_filename, 'w', driver='GTiff', width=1000, height=1000, count=1,
                           dtype='uint8') as dst:
            dst.write(image, 1)
            dst.write_colormap(

           1, {

        0: (229, 245, 249),
        1: (214, 238, 236),
        2: (198, 231, 223),
        3: (183, 224, 211),
        4: (167, 217, 198),
        5: (152, 210, 185),
        6: (136, 203, 172),
        7: (121, 197, 159),
        8: (106, 190, 146),
        9: (90, 183, 133),
        10: (75, 176, 121),
        11: (59, 169, 108),
        12: (44, 162, 95),
        13: (44, 162, 95),
        14: (44, 162, 95),
        15: (44, 162, 95),
        16: (44, 162, 95),
        17: (44, 162, 95),
        18: (44, 162, 95),
        19: (44, 162, 95),
        20: (44, 162, 95),
        21: (44, 162, 95),
        22: (44, 162, 95),
        23: (44, 162, 95),
        24: (44, 162, 95),
        25: (44, 162, 95) })
        self.temp_files.append(out_filename)

    def to_video(self,path):
        result_code = subprocess.call("ffmpeg -framerate 25 -pattern_type glob -i '{}*_state.tif'  {}output.mp4 -y".format(path,path), shell=True)

    def clean_up(self):
        while self.temp_files:
            file = self.temp_files.pop()
            os.remove(file)
            print("removed tmp file {}".format(file))
            
if __name__ == "__main__":
    import glob
    import os

    root = 'output'
    plotter = Plotter()

    for p in os.listdir(root):
        path = os.path.join(root,p+"/")
        print("searching in {}".format(path))
        paths = glob.glob(path+'*_state.csv')
        print('paths = ',paths)
        for image_path in paths:
            print('path = ',image_path)
            plotter.reconstruct(image_path)
 #       plotter.to_video(path)
 #       plotter.clean_up()




