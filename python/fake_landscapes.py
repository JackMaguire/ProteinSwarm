#import numpy as np
from random import *
import math
import numpy as np

def add_noise( score, noise_level=0.01 ):
    one_percent = score * noise_level
    noise = uniform( 0, one_percent )
    return score * ( 1 + noise )

class FakeLandscape:
    #This landscape will have multiple minima of random depths

    def __init__(self, ndim=6, noise_level=0.05, width_factor = 1.0):
        self.ndim = ndim #Expecting 1 to 6 dimensions

        #min and max values for minima
        self.min = -1.0
        self.max = 1.0
        self.scale = 1.0 #we want to stay between 0 and 100

        #define min and max depths for the wells
        #(max means more negative, I suppose)
        self.max_depth = -1
        self.min_depth = -0.5

        self.noise = noise_level
        
        self.n_minima = 3 * self.ndim
        self.minima = [] #each element is an array
        self.minima_depths = []
        self.minima_factor = [] #make the wells more shallow
        self.global_minimum = 0

        for i in range( 0, self.n_minima ):
            minima_coords = []
            for _ in range( 0, self.ndim ):
                minima_coords.append( uniform( self.min, self.max ) )

            self.minima.append( minima_coords )
            self.minima_factor.append( uniform( 0.01, 0.1 ) * width_factor )
            depth = uniform( self.max_depth, self.min_depth )
            self.minima_depths.append( depth )
            if depth < self.global_minimum:
                self.global_minimum = depth


    #arr is a numpy array with shape (self.ndim)
    def score( self, arr, noise=True ):
        #scale
        #arr := {-1,1}
        #landscape := {-100,100}
        vals = []
        for a in arr:
            vals.append( a * self.scale )

        #let's just do a simple quadratic to the closest well
        #Energy wells are negative, let's say the lanscape is flat at +1
        best_score = 1
        for i in range( 0, self.n_minima ):
            distance = 0
            for j in range( 0, self.ndim ):
                distance += (vals[ j ] - self.minima[ i ][ j ]) ** 2
            distance = math.sqrt( distance )

            quardtratic = ( distance * self.minima_factor[ i ] ) ** 2
            score = quardtratic + self.minima_depths[ i ]
            if score < best_score:
                best_score = score
        if noise:
            return add_noise( best_score, self.noise )
        else:
            return best_score
