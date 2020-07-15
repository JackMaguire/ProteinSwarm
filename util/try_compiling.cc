#include "../include/protein_swarm.hh"
#include <vector>

using namespace protein_swarm;

int main(){
  std::vector< Bounds > bounds( 2 );

  bounds[ 0 ].lower_bound = -1;
  bounds[ 0 ].upper_bound = 1;
  bounds[ 0 ].type = BoundsType::STANDARD;

  bounds[ 1 ].lower_bound = 0;
  bounds[ 1 ].upper_bound = 10;
  bounds[ 1 ].type = BoundsType::PACMAN;

  protein_swarm::ProteinSwarm swarm ( 5, 2, bounds );
  auto x = swarm.ask( 0.0 );
  swarm.tell( x, -1.0 );
}
