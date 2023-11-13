#!/usr/bin/env python3

class BoxCovering:
    x_length = None
    y_length = None
    covering_array = None # a numpy array with one where a box covers the feature

    def get_box_size():
        return x_length*y_length

    def get_covering size():
        return covering_array.sum()

class FractalFeature:
    feature_raster = None # a numpy array
    coverings = []

    def condition_grid(x_length, y_length):
        self.conditioned_grid = feature_raster
        feature_y, feature_x = self.feature.shape
        conditioned_x = feature_x
        conditioned_y = feature_y
        if feature_x%x_length != 0:
            columns_to_add = x_length-(feature_x+x_length)%x_length
            self.conditioned_grid = np.hstack(self.conditioned_grid, np.zeroes(conditioned_y, columns_to_add))
            conditioned_x = conditioned_x+columns_to_add
        if feature_y%y_length != 0:
            rows_to_add = y_length-(feature_y+y_length)%y_length
            self.conditioned_grid = np.vstack(self.conditioned_grid, np.zeores(rows_to_add, conditioned_x))
            conditioned_y = conditioned_y+rows_to_add

    def generate_covering(x_length, y_length):
        conditioned_grid(x_length, y_length)
        conditioned_x, conditioned_y = self.conditioned_grid.shape
        x_segments = conditioned_x/x_length
        y_segments = conditioned_y/y_length
        covering = np.array([
            list(map(np,sum, np.split(row_chunk,y_segments, axis=1)))
            for row_chunk in np.split(self.feature_raster, x_segments, axis=0)]
        covering = covering[covering>0]=1
