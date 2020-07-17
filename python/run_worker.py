from mpi4py import MPI
from fake_landscapes import *
#import nevergrad as ng
#from numpy import mean
import numpy as np
import random as rand
#import argparse

from proteinswarm import *


def score( dofs ):
    assert( dofs.shape[0] == 7 )
    total_score = 0.0
    seed = 0
    for ndim in [ 10, 15, 20 ]:
        for width_factor in [ 4, 8, 16, 32, 64 ]:
            for budget_factor in [ 500, 750, 1000 ]:
                seed += 1
                np.random.seed( seed )
                rand.seed( a=seed )

                budget=ndim*budget_factor

                landscape = FakeLandscape( ndim=ndim, noise_level=0.05, width_factor=width_factor )
                bounds = ProteinSwarmBounds( ndim )
                bounds.set_all_lower_bounds( -3 )
                bounds.set_all_upper_bounds( 3 )
                optimizer = bounds.make_pso( int(budget_factor/3) )

                optimizer.set_span_coeff( dofs[ 0 ] )
                optimizer.set_k1( dofs[ 1 ] )
                optimizer.set_k2( dofs[ 2 ] )
                optimizer.set_c1( dofs[ 3 ] )
                optimizer.set_c2( dofs[ 4 ] )
                optimizer.set_v_limit_m( dofs[ 5 ] )
                optimizer.set_v_limit_b( dofs[ 6 ] )

                if not optimizer.parameters_are_reasonable():
                    return 10.0

                for b in range( 0, budget ):
                    sample = optimizer.ask( float(b) / float(budget) )
                    score = landscape.score( extract_data( sample ) )
                    optimizer.tell( sample, score )

                best_score = optimizer.get_global_best_score() / landscape.global_minimum
                total_score += best_score

    return -1.0 * (total_score / float( seed ))

def run_worker( comm, rank ):

    while True:
        status = MPI.Status()
        dofs = comm.recv( source=0, tag=MPI.ANY_TAG, status=status )
        if status.Get_tag() == 0:
            comm.send( 0, dest=0, tag=0 )
            break

        final_score = score( dofs.value )
        bundle = [ dofs, final_score ]
        comm.send( bundle, dest=0, tag=1 )
