from mpi4py import MPI
from fake_landscape import *
import nevergrad as ng
#from numpy import mean
import numpy as np
import random as rand
#import argparse

import proteinswarm as pso

def score( dofs ):
    total_score = 0
    seed = 0
    for ndim in [ 10, 15, 20, 25, 30 ]:
        for width_factor in [ 4, 8, 16, 32, 64 ]:
            for noise in [ 0.0, 0.05, 0.1, 0.15 ]:
                seed += 1
                np.random.seed( seed )
                rand.seed( a=seed )

                landscape = FakeLandscape( ndim=ndim, noise_level=noise, width_factor=width_factor )
                # TODO


    return total_score

def run_worker( comm, rank ):

    while True:
        status = MPI.Status()
        dofs = comm.recv( source=0, tag=MPI.ANY_TAG, status=status )
        if status.Get_tag() == 0:
            comm.send( 0, dest=0, tag=0 )
            break

        final_score = score_dofs( dofs )
        bundle = [ dofs, final_score ]
        comm.send( bundle, dest=0, tag=1 )
