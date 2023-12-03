#!/usr/bin/env python3
import numpy as np
from osgeo import gdal
import os
from pathlib import Path

def covering_score(covering, x_length, y_length):
        """
        This "scores" a covering based on how close below it is the size of a new deisred covering
        """
        covering_score = covering.x_length/x_length+covering.y_length/y_length
        if covering_score>2:
            covering_score = -1
        return covering_score

class BoxCovering:
    def __init__(self, covering_array, x_length, y_length):
        self.x_length = None
        self.y_length = None
        self.covering_array = None # a numpy array with one where a box covers the feature

    def get_box_size(self):
        """
        Returns the size of the boxes that are used in the covering.  This can be in any units
        but is typically in pixels of the initial feature raster.
        """
        return self.x_length*self.y_length

    def get_covering_size(self):
        """
        Returns the total number of boxes in the covering
        """
        return np.sum(self.covering_array)

class GdalFractal:
    def __init__(self, shapefile_path, raster_directory, initial_res):
        self.shapefile_path = shapefile_path
        self.raster_directory = raster_directory
        Path(raster_directory).mkdir(parents=True, exist_ok=True)
        self.raster_format='GTIFF'
        self.creation_options=["COMPRESS=DEFLATE"]
        self.covering_points = []
        self.add_covering(initial_res)

    def add_covering(self, resolution):
        covering_name = "fractal_%d.tif" % resolution
        covering_path = os.path.join(self.raster_directory, covering_name)
        covering_object = gdal.Rasterize(covering_path, self.shapefile_path, format=self.raster_format,
                                  creationOptions=self.creation_options, noData=0, initValues=1,
                                  xRes=resolution, yRes=-resolution, allTouched=True, burnValues=1)
        covering_object = None
        print(covering_path)
        covering = gdal.Open(covering_path)
        covering_array = covering.ReadAsArray()
        covering = None
        self.covering_points.append((resolution, np.sum(covering_array)))
        covering_array = None

    def return_covering_points(self):
        box_lengths = np.array([p[0] for p in self.covering_points])
        box_numbers = np.array([p[1] for p in self.covering_points])
        return box_lengths, box_numbers


class FractalFeature:
    def __init__(self, feature, x_length=1, y_length=1):
        self.feature_raster = feature# a numpy array
        self.x_length=x_length
        self.y_length=y_length
        self.coverings = [self]
        self.covering_points = [(np.sqrt(x_length*y_length), np.sum(feature))]

    def return_covering_points(self):
        box_lengths = np.array([p[0] for p in self.covering_points])
        box_numbers = np.array([p[1] for p in self.covering_points])
        return box_lengths, box_numbers

    def condition_grid(self, x_length, y_length):
        """
        This adds rows and columns of zeroes to make the grid is evenly visibile by the size of the covering boxes
        """
        #x_length= x_length//self.x_length
        #y_length = y_length//self.y_length
        self.conditioned_grid = self.feature_raster
        feature_y, feature_x = self.feature_raster.shape
        conditioned_x = feature_x
        conditioned_y = feature_y
        if feature_x%x_length != 0:
            columns_to_add = x_length-(feature_x+x_length)%x_length
            self.conditioned_grid = np.hstack((self.conditioned_grid, np.zeros((conditioned_y, columns_to_add))))
            conditioned_x = conditioned_x+columns_to_add
        if feature_y%y_length != 0:
            rows_to_add = y_length-(feature_y+y_length)%y_length
            self.conditioned_grid = np.vstack((self.conditioned_grid, np.zeros((rows_to_add, conditioned_x))))
            conditioned_y = conditioned_y+rows_to_add

    def covering_score(covering, x_length, y_length):
        """
        This "scores" a covering based on how close below it is the size of a new deisred covering
        """
        covering_score = covering.x_length/x_length+covering.y_length/y_length
        if covering_score>2:
            covering_score = -1
        return covering_score

    def get_closest_covering(self, x_length, y_length):
        """
        Returns the cloest below covering to the size of a desired new covering
        """
        closest_covering = max(self.coverings, key=lambda c: covering_score(c, x_length, y_length))
        return closest_covering

    def get_valid_subcoverings(self, x_length, y_length):
        """
        Returns only coverings of sizes that are proper divisors of the desired covering size
        """
        return [covering for covering in self.coverings if ((x_length%covering.x_length) + (y_length%covering.y_length))==0]

    def get_closest_valid_subcovering(self, x_length, y_length):
        """
        Returns the closest below covering that is a proper divisor of the deisred covering size
        """
        closest_covering = max(self.get_valid_subcoverings(x_length, y_length), key=lambda c: covering_score(c, x_length, y_length))
        return closest_covering

    def generate_covering(self, x_length, y_length):
        """
        Generates a new fractal object of the desired size
        """
        x_length= x_length//self.x_length
        y_length = y_length//self.y_length
        self.condition_grid(x_length, y_length)
        conditioned_y, conditioned_x = self.conditioned_grid.shape
        x_segments = conditioned_x/x_length
        y_segments = conditioned_y/y_length
        covering = np.array([
            list(map(np.sum, np.split(row_chunk,x_segments, axis=1)))
            for row_chunk in np.split(self.conditioned_grid, y_segments, axis=0)])
        covering[covering>0]=1
        covering_object= FractalFeature(covering, x_length, y_length)
        return covering_object

    def add_covering(self, covering=None, x_length=None, y_length=None, existing=True):
        """ adds a covering to the list of coverings.  Either a provided covering or a new one is generated"""
        if covering is not None:
            self.coverings.append(covering)
        elif x_length is not None:
            if y_length is None:
                y_length=x_length
            if existing:
                covering = self.generate_covering_from_existing(x_length, y_length)
            else:
                covering = self.generate_covering(x_length, y_length)
            self.coverings.append(covering)

        else:
            raise Exception("No covering or generating parameters provided")
        self.covering_points.append((np.sqrt(covering.x_length*covering.y_length), np.sum(covering.feature_raster)))

    def generate_covering_from_existing(self, x_length, y_length):
        """
        Generates a new fractal object of the desired size from the closest existing covering
        """
        closest_covering = self.get_closest_valid_subcovering(x_length, y_length)
        if closest_covering==self:
            covering = self.generate_covering(x_length, y_length)
        else:
            #x_length= x_length/closest_covering.x_length
            #y_length = y_length/closest_covering.y_length
            covering = closest_covering.generate_covering_from_existing(x_length, y_length)
        return covering
