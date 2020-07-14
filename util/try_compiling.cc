#include "../include/protein_swarm.hh"

int main(){
  protein_swarm::ProteinSwarm< float > swarm ( 5, 2, {0,1}, {-1,1} );
  auto x = swarm.ask( 0.0 );
  swarm.tell( x, -1.0 );
}
