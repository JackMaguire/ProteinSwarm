//g++ unit_tests.cc -o unit_tests -Wall -std=c++11
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


int main(){
  run_standard_bounds_tests();
  run_pacman_bounds_tests();
  std::cout << "All tests passed!" << std::endl;
}
