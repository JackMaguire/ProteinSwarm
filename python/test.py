from proteinswarm import *

from fake_landscapes import *

ndim = 20

bounds = ProteinSwarmBounds( ndim )
bounds.set_all_lower_bounds( -1 )
bounds.set_all_upper_bounds( 1 )
bounds.set_all_bounds_type( "standard" )
#bounds.set_bounds_type( 0, "pacman" )

nparticles = 10
optimizer = bounds.make_pso( nparticles )

landscape = FakeLandscape( ndim=ndim, noise_level=0.1, width_factor = 30.0 )

print( "done setting up" )

for x in range( 0, 1000 ):
    sample = optimizer.ask( float(x) / 1000.0 )
    score = landscape.score( extract_data( sample ) )
    print( x, score )
    optimizer.tell( sample, score )
