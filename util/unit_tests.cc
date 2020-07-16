//g++ unit_tests.cc -o unit_tests -Wall -std=c++11 -D_GLIBCXX_DEBUG
#define UNIT_TEST

#include "../include/protein_swarm.hh"

//This second include is a test in itself
#include "../include/protein_swarm.hh"

#include <vector>
#include <assert.h>
#include <iostream>

using namespace protein_swarm;

void run_standard_bounds_tests(){
  std::cout << " - running run_standard_bounds_tests()" << std::endl;
  {
    Bounds b;
    b.lower_bound = -1.0;
    b.upper_bound = 1.0;
    b.type = BoundsType::STANDARD;

    b.verify();

    assert( b.span() == 2.0 );

    assert( b.apply( -2.0 ) == -1.0 );
    assert( b.apply( -1.0 ) == -1.0 );
    assert( b.apply( 0.0 ) == 0.0 );
    assert( b.apply( 1.0 ) == 1.0 );
    assert( b.apply( 2.0 ) == 1.0 );

    //inside range
    assert( b.get_vector( -0.5, 0.5 ) == 1.0 );
    assert( b.get_vector( 0.25, -0.5 ) == -0.75 );
    assert( b.get_vector( 0.25, 0.25 ) == 0 );

    //outside range - this should still work
    assert( b.get_vector( -100.0, 0 ) == 100.0 );
    assert( b.get_vector( 500.0, -2.0 ) == -502.0 );
  }

  {
    Bounds b;
    b.lower_bound = 0.0;
    b.upper_bound = 10.5;
    b.type = BoundsType::STANDARD;

    b.verify();

    assert( b.span() == 10.5 );

    assert( b.apply( -1.0 ) == 0.0 );
    assert( b.apply(  0.0 ) == 0.0 );
    assert( b.apply( 10.0 ) == 10.0 );
    assert( b.apply( 10.5 ) == 10.5 );
    assert( b.apply( 11.0 ) == 10.5 );
  }

  //Negative Tests
  {
    Bounds b;
    b.lower_bound = 0;
    b.upper_bound = 0;
    //This should fail because we need lower_bound < upper_bound

    bool threw = false;
    try {
      b.verify();
    } catch(...) {
      threw = true;
    }
    assert( threw );
  }

  {
    Bounds b;
    b.lower_bound = -1.0;
    b.upper_bound = -2.0;
    //This should fail because we need lower_bound < upper_bound

    bool threw = false;
    try {
      b.verify();
    } catch(...) {
      threw = true;
    }
    assert( threw );
  }

}

void run_pacman_bounds_tests(){
  std::cout << " - running run_pacman_bounds_tests()" << std::endl;
  {
    Bounds b;
    b.lower_bound = -1.0;
    b.upper_bound = 1.0;
    b.type = BoundsType::PACMAN;

    b.verify();

    assert( b.span() == 2.0 );

    assert( b.apply( -2.0 ) == 0.0 );
    assert( b.apply( -1.5 ) == 0.5 );
    assert( b.apply( -1.0 ) == -1.0 );
    assert( b.apply( 0.0 ) == 0.0 );
    assert( b.apply( 1.0 ) == 1.0 );
    assert( b.apply( 1.25 ) == -0.75 );
    assert( b.apply( 2.0 ) == 0.0 );

    //inside range
    assert( b.get_vector( 0, 1.5 ) == -0.5 );
    assert( b.get_vector( 0.25, -0.5 ) == -0.75 );
    assert( b.get_vector( -0.25, 0.5 ) == 0.75 );
    assert( b.get_vector( 0.75, -0.5 ) == 0.75 );
    assert( b.get_vector( 0.75, 0.75 ) == 0.0 );
    assert( b.get_vector( 0.5, 1.0 ) == 0.5 );
    assert( b.get_vector( 0.5, -1.0 ) == 0.5 );

    //outside range - this should still work?
    assert( b.get_vector( -100.0, -100.0 ) == 0.0 );    
    assert( b.get_vector( -1.25, -1.0 ) == 0.25 );
    assert( b.get_vector( -1.25, 1.0 ) == 0.25 );
    assert( b.get_vector( 1.25, 0.0 ) == 0.75 );
    assert( b.get_vector( 1.25, 1.0 ) == -0.25 );

  }

  {
    Bounds b;
    b.lower_bound = 0.0;
    b.upper_bound = 10.5;
    b.type = BoundsType::PACMAN;

    b.verify();

    assert( b.span() == 10.5 );

    assert( b.apply( -1.0 ) == 9.5 );
    assert( b.apply(  0.0 ) == 0.0 );
    assert( b.apply( 10.0 ) == 10.0 );
    assert( b.apply( 10.5 ) == 10.5 );
    assert( b.apply( 11.0 ) == 0.5 );
  }

}

void run_simple_ask_tell(){
  std::cout << " - running run_simple_ask_tell()" << std::endl;  

  std::vector< Bounds > bounds( 1 );
  bounds[ 0 ].lower_bound = 0.0;
  bounds[ 0 ].upper_bound = 1.0;

  ProteinSwarm optimizer( 2, bounds );

  assert( optimizer.get_n_particles_in_queue() == 2 );

  Sample const x1 = optimizer.ask( 0.0 );
  assert( x1.info.particle == 0 );

  assert( optimizer.get_n_particles_in_queue() == 1 );

  Sample const x2 = optimizer.ask( 0.0 );
  assert( x2.info.particle == 1 );

  assert( optimizer.get_n_particles_in_queue() == 0 );

  optimizer.tell( x1, -1.0 );//NEW BEST
  assert( optimizer.get_index_of_global_best() == 0 );
  assert( optimizer.get_global_best_position().size() == 1 );
  assert( optimizer.get_global_best_position()[ 0 ] == x1.value[ 0 ] );

  assert( optimizer.get_n_particles_in_queue() == 1 );

  optimizer.tell( x2, -2.0 );//NEW BEST
  assert( optimizer.get_index_of_global_best() == 1 );
  assert( optimizer.get_global_best_position().size() == 1 );
  assert( optimizer.get_global_best_position()[ 0 ] == x2.value[ 0 ] );

  assert( optimizer.get_n_particles_in_queue() == 2 );

  Sample const x3 = optimizer.ask( 0.0 );
  assert( x3.info.particle == 0 );
  optimizer.tell( x2, 0.0 );
  assert( optimizer.get_index_of_global_best() == 1 ); //NO CHANGE
  assert( optimizer.get_global_best_position()[ 0 ] == x2.value[ 0 ] ); //NO CHANGE

  assert( optimizer.get_n_particles_in_queue() == 2 );
}

int main(){
  run_standard_bounds_tests();
  run_pacman_bounds_tests();
  run_simple_ask_tell();
  std::cout << "All tests passed!" << std::endl;
}
