#!/usr/bin/env python3

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

class FractalFeature:
    def __init__(self, feature, x_length=1, y_length=1):
        feature_raster = None # a numpy array
        self.x_length=None
        self.y_length=None
        coverings = [self]

    def condition_grid(self, x_length, y_length):
        """
        This adds rows and columns of zeroes to make the grid is evenly visibile by the size of the covering boxes
        """
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
        return [covering for covering in self.coverings if (covering.x_length%x_length) + (covering.y_length%y_length)==0]

    def get_closest_valid_subcovering(self):
        """
        Returns the closest below covering that is a proper divisor of the deisred covering size
        """
        closest_covering = max(get_valid_subcoverings(self, key=lambda c: covering_score(c, x_length, y_length)))
        return closest_covering

    def generate_covering(self, x_length, y_length):
        """
        Generates a new fractal object of the desired size
        """
        condition_grid(x_length, y_length)
        conditioned_x, conditioned_y = self.conditioned_grid.shape
        x_segments = conditioned_x/x_length
        y_segments = conditioned_y/y_length
        covering = np.array([
            list(map(np,sum, np.split(row_chunk,y_segments, axis=1)))
            for row_chunk in np.split(self.feature_raster, x_segments, axis=0)])
        covering = covering[covering>0]=1
        return covering

    def generate_covering_from_existing(self, x_length, y_length):
        """
        Generates a new fractal object of the desired size from the closest existing covering
        """
        closest_covering = get_closest_valid_covering(x_length, y_length)
        if cloest_covering==self:
            covering = generate_covering(x_length, y_length)
        else:
            x_length= x_length/closest_covering.x_length
            y_length = y_length/closest_covering.y_length
            covering = closest_covering.generat_covering_from_existing(x_length, y_length)
        self.coverings.append(covering)
