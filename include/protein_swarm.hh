// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-

#pragma once

#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <assert.h>
#include <limits>
#include <cstdlib> //rand()
#include <math.h>
#include <chrono>

//#include <array>
#include <vector>
#include <queue>

/*
  CURRENT ASSUMPTIONS
  - N Partcles >= N CPUs. There will never be more than one CPU working on the same particle at the same time
	- The cost of the simulation >> the cost of running this code. This is not optimized for runtime performance.
	- You are currently required to provide bounds
	- This is not internally thread safe, you may want to wrap a mutex around your ProteinSwarm object
*/

namespace protein_swarm {

using uint = unsigned int;
using ulong = unsigned long int;
using Value = double;


namespace {

Value rand_01(){
	return Value( rand() ) / Value( RAND_MAX );
}

#ifdef UNIT_TEST
//Unit tests will sometimes need to assert that an assertion is thrown

//This hack gets the job done

struct ProteinSwarmError{};

void unit_test_assert( bool b ){
	if( ! b ) throw ProteinSwarmError{};
}

#undef assert

#define assert( b ) unit_test_assert( b )
#endif


} //anonymous namespace



enum class BoundsType
	{
	 STANDARD, //Can't go above the upper bound or below the lower bound

	 PACMAN //Once you cross the boundary, you appear on the other side.
	};



struct Bounds {
	Value lower_bound;
	Value upper_bound;
	BoundsType type = BoundsType::STANDARD;

	Value span() const {
		return upper_bound - lower_bound;
	}

	Value get_vector( Value source, Value destination ) const;

	Value apply( Value unbounded ) const;

	//assert that the current values are well behaved
	void verify() const {
		assert( ! isnan( lower_bound ) );
		assert( ! isnan( upper_bound ) );
		assert( lower_bound < upper_bound );
	};
};



//Commenting out enums that are not implemented yet
enum class InitialSamplingMethod
  {

   //every particle randomly chooses its starting position from [ lower_bound, upper_bound ] uniformly
   UNIFORM,

   //every particle randomly chooses its starting position from [ -1, 1 ] uniformly
   //UNIFORM_NORM,

   //every particle randomly chooses its starting position from [ 0, 1 ] uniformly
   //UNIFORM_NORM_0,

   //every particle randomly chooses its starting position by sampling from a gaussian distribution centered at lower_bound + (upper_bound-lower_bound)/2
   //with a std deviation of (upper_bound-lower_bound)/6 and clipped at [ lower_bound, upper_bound ]
   //GAUSSIAN,

   //every particle randomly chooses its starting position by sampling from a gaussian distribution centered at 0
   //with a std deviation of (upper_bound-lower_bound)/6 and clipped at [ lower_bound, upper_bound ]
   //GAUSSIAN_0,

   //every particle randomly chooses its starting position by sampling from a gaussian distribution centered at 0
   //with a std deviation of 1 and clipped at [ lower_bound, upper_bound ]
   //GAUSSIAN_NORM,

   //particles attempt to be evenly spaced out
   //GRID,

   //use user-provided intialization by calling set_particle_initializer
   //CUSTOM
  };



class ParticleInitialzer {
public:
  virtual ~ParticleInitialzer() = default;

  //This will be called many times for a single run, once for each particle
  virtual
  std::vector< Value >
  initialize_one_particle(
		uint ndim,//How many dimensions does this case have?
    uint n_total_particles, //How many particles are in the system?
    uint particle_id, //which particle are we considering now? Ranging from [ 0, n_total_particles-1 ] inclusive
    std::vector< Bounds > const & bounds//these are the bounds, which the user probably already knows anyways
  ) const = 0; //returns the starting values for each dimension for this particle

};



class UniformParticleInitialzer : public ParticleInitialzer {
public:
	UniformParticleInitialzer() = default;

  ~UniformParticleInitialzer() override = default;

  //This will be called many times for a single run, once for each particle
  std::vector< Value >
  initialize_one_particle(
		uint ndim,
    uint,
    uint,
    std::vector< Bounds > const & bounds
  ) const override;
};

//This is one case where array of struct appears to be better than struct of array
class Particle {

public:
	Particle( Particle const & ) = default;
	Particle( Particle && ) = default;

	Particle(){
		//current_position_.resize( ndim ); //keep positions empty until they are set
		//current_speed_.resize( ndim ); //keep empty until it is set
		//best_position_.resize( ndim ); //keep positions empty until they are set
	}

	~Particle() = default;

public:
	void
	initialize(
		std::vector< Value > const & current_position,
		std::vector< Value > const & current_speed
	){
		current_position_ = current_position;
		current_speed_ = current_speed;
		best_position_ = current_position;
	}

public: //accessors
	void set_current_score_is_outdated( bool setting ){
		current_score_is_outdated_ = setting;
	}

	Value get_current_score() const {
		assert( ! current_score_is_outdated_ );
		return current_score_;
	}

	void set_current_score( Value v ){
		current_score_ = v;
		current_score_is_outdated_ = false;
		if( current_score_ < best_score_ ){
			best_position_ = current_position_;
		}
	}

	std::vector< Value > const & get_current_position() const {
		return current_position_;
	}

	void set_current_position( std::vector< Value > const & v ){
		current_position_ = v;
	}

	void set_current_position( std::vector< Value > && v ){
		current_position_ = v;
	}

	void set_current_position( uint const dim, Value const setting ){
		current_position_[ dim ] = setting;
	}

	std::vector< Value > const & get_current_speed() const {
		return current_speed_;
	}

	void set_current_speed( std::vector< Value > const & v ){
		current_speed_ = v;
	}

	void set_current_speed( std::vector< Value > && v ){
		current_speed_ = v;
	}

	void set_current_speed( uint const dim, Value const setting ){
		current_speed_[ dim ] = setting;
	}

	Value get_best_score() const {
		return best_score_;
	}

	void set_best_score( Value v ){
		best_score_ = v;
	}

	std::vector< Value > const & get_best_position() const {
		return best_position_;
	}

	void set_best_position( std::vector< Value > const & v ){
		best_position_ = v;
	}

	void set_best_position( std::vector< Value > && v ){
		best_position_ = v;
	}

private:
	Value current_score_ = std::numeric_limits< Value >::max();
	std::vector< Value > current_position_; //Keep empty until it is set
	std::vector< Value > current_speed_;

	Value best_score_ = std::numeric_limits< Value >::max();
	std::vector< Value > best_position_; //Keep empty until it is set

	bool current_score_is_outdated_ = true;
	//This is true when the current_position_ has been set to the new value but the current_score_ has not been calculated yet. current_score_is_outdated_ == true means that the job is currently in flight

};

struct SampleInfo {
  uint particle;
};

struct Sample {
	SampleInfo info;
	std::vector< Value > value;
};

//THIS IS NOT THREADSAFE
class ProteinSwarm {
	static_assert( std::is_floating_point< Value >::value, "Protein swarm was written assuming floating point values." );

public:

  ProteinSwarm(
    uint n_particles,
		//uint ndim,
    std::vector< Bounds > const & bounds,
    InitialSamplingMethod sampling_method = InitialSamplingMethod::UNIFORM
  );

  void reset(){
    //TODO
  }

	Sample ask( float fraction_of_run_completed );

  void tell( Sample const & sample, Value score ){
    tell( sample.info, score );
  }

  void tell( SampleInfo info, Value score );

public://getters and setters
  uint get_n_particles() const {
    return n_particles_;
  }

  /*void set_n_particles( uint setting ){
    n_particles_ = setting;
  }*/

	Value get_global_best_score() const {
		return global_best_score_;
	}

	std::vector< Value > const &
	get_global_best_position() const {
		return particles_[ index_of_global_best_ ].get_best_position();
	}

public://unit test
	uint get_index_of_global_best() const {
		return index_of_global_best_;
	}

	uint get_n_particles_in_queue() const {
		return particle_queue_.size();
	}

protected:
  void initialize( ParticleInitialzer const & initializer );

	void
	update_to_new_position(
		uint particle_id,
		float fraction_of_run_completed // [0.0, 1.0)
	);

private:
  uint n_particles_;
	uint ndim_;

	//Particle holds location, history, and score info for each particle
	std::vector< Particle > particles_;

  uint index_of_global_best_ = 0;
	Value global_best_score_ = std::numeric_limits< Value >::max();

	//The particle queue holds a list of particles that are ready to be sampled
	std::queue< uint > particle_queue_;

	std::vector< Bounds > bounds_;

  uint n_asks_ = 0;
	uint n_tells_ = 0;
};



inline
Value
Bounds::get_vector( Value source, Value destination ) const {
	switch( type ){
	case( BoundsType::STANDARD ):
		return destination - source;
	case( BoundsType::PACMAN ):
		//Example 1:
		//lower_bound = -1
		//upper_bound = 1
		//source = 0.5
		//destination = -0.75
		//dist = 1.25
		//pacman_dist = 0.75
		//vector = 0.75

		//Example 2:
		//lower_bound = -1
		//upper_bound = 1
		//source = -0.5
		//destination = 0.75
		//dist = 1.25
		//pacman_dist = 0.75
		//vector = -0.75

		//Apply transformations just in case
		source = apply( source );
		destination = apply( destination );

		Value const dist = abs( destination - source );
		Value const pacman_dist = span() - dist;

		if( dist < pacman_dist ) return destination - source;

		//else, use pacman
		if( destination < source ) //then we want to go in the positive direction
			return pacman_dist;
		else //negative direction
			return -pacman_dist;
	}
	assert( false ); //dead code
}



  // // // //
// FUNCTIONS //
  // // // //

inline
Value
Bounds::apply( Value unbounded ) const {
	switch( type ){
	case( BoundsType::STANDARD ):
		return std::min( std::max( unbounded, lower_bound ),	upper_bound );

	case( BoundsType::PACMAN ):
		Value const s = span();
		//using while instead of modulo assuming that we are usually not too far out of bounds
		while( unbounded > upper_bound ) unbounded -= s;
		while( unbounded < lower_bound ) unbounded += s;
		return unbounded;
	}
	assert( false ); //dead code
}

inline
std::vector< Value >
UniformParticleInitialzer::initialize_one_particle(
	uint ndim, uint, uint,
	std::vector< Bounds > const & bounds
) const {
	std::vector< Value > vec( ndim );
	for( uint d = 0; d < ndim; ++d ){
		vec[ d ] = bounds[ d ].lower_bound + ( bounds[d].span() * rand_01() );
	}
	return vec;
}

inline
ProteinSwarm::ProteinSwarm(
	uint const n_particles,
	//uint const ndim,
	std::vector< Bounds > const & bounds,
	InitialSamplingMethod const sampling_method
):
	n_particles_( n_particles ),
	ndim_( bounds.size() ),
	bounds_( bounds )
{
	unsigned long microseconds_since_epoch =
    std::chrono::duration_cast< std::chrono::microseconds >
		(std::chrono::system_clock::now().time_since_epoch()).count();
	srand( microseconds_since_epoch );
	assert( bounds_.size() == ndim_ );

	switch( sampling_method ){
	case( InitialSamplingMethod::UNIFORM ):
		UniformParticleInitialzer init; //TODO save as unique_ptr?
		initialize( init );
		break;
	}
}

inline
Sample
ProteinSwarm::ask( float const fraction_of_run_completed ){
	++n_asks_;

	assert( ! particle_queue_.empty() );
	uint const particle = particle_queue_.front();
	particle_queue_.pop();

	//if this is true then every particle has been scored at least once
	if( n_asks_ > n_particles_ ){
		//Move the particle if this is not its initial measurement
		update_to_new_position(	particle,	fraction_of_run_completed );
	}

	Sample s;
	s.info.particle = particle;
	s.value = particles_[ particle ].get_current_position();
	return s;
}

inline
void
ProteinSwarm::tell( SampleInfo const info, Value const score ){
	++n_tells_;
	particles_[ info.particle ].set_current_score( score );

	if( score < global_best_score_ ){
		index_of_global_best_ = info.particle;
		global_best_score_ = score;
	}

	particle_queue_.push( info.particle );
}

inline
void
ProteinSwarm::initialize(
	ParticleInitialzer const & initializer
){
	assert( ndim_ > 0 );
	assert( n_particles_ > 0 );
	assert( ! bounds_.empty() );

	for( Bounds const & b : bounds_ ){
		b.verify();
	}

	particles_.resize( n_particles_ );
	for( uint i = 0; i < n_particles_; ++i ){
		particle_queue_.push( i );

		std::vector< Value > const starting_position =
			initializer.initialize_one_particle( ndim_, n_particles_, i, bounds_ );

		std::vector< Value > starting_velocity( ndim_ );
		for( uint d = 0; d < ndim_; ++d ){
			constexpr Value span_coeff = 0.25;
			Value const vec = ( bounds_[ d ].span() * span_coeff ) //scale magnitude of V to range size
				* ( 2*rand_01() - 1.0); // adjust to ( -25% to 25% )
			starting_velocity[ d ] = vec;
		}

		particles_[ i ].initialize( starting_position, starting_velocity );
	}//for i

}//initialize

inline
void
ProteinSwarm::update_to_new_position(
	uint const particle_id,
	float const fraction_of_run_completed // [0.0, 1.0)
){
	Particle & particle = particles_[ particle_id ];

	assert( ! particle.get_current_position().empty() ); //Did you call initialize()?
	assert( particle.get_current_speed().size()  == particle.get_current_position().size() );
	assert( particle.get_best_position().size() == particle.get_current_position().size() );


	std::vector< Value > const & global_best_position = particles_[ index_of_global_best_ ].get_best_position();

	assert( global_best_position.size() == particle.get_current_position().size() );
	assert( bounds_.size() == particle.get_current_position().size() );

	particle.set_current_score_is_outdated( true );

	static_assert( std::is_floating_point< Value >::value, "Protein swarm was written assuming floating point values." );

	//w = 0.9 - ( (0.7 * t) / numofiterations);
	//w = k1 - ( k2 * fraction_of_run_completed);
	constexpr Value k1 = 0.9;
	constexpr Value k2 = 0.7;
	Value const W = k1 - ( k2 * fraction_of_run_completed );

	constexpr Value c1 = 2.0; //not sure what to think of this
	constexpr Value c2 = 2.0; //not sure what to think of this

	//sample a new position
	for( uint d = 0; d < ndim_; ++d ){
		Value const curr_pos_d = particle.get_current_position()[ d ];
		Value const best_pos_d = particle.get_best_position()[ d ];

		//V[i][j] = min(max((w * V[i][j] + rand_01 * c1 * (pbests[i][j] - X[i][j]) + rand_01 * c2 * (gbest[j] - X[i][j])), Vmin[j]), Vmax[j]);
		Value const vec_to_local_best =
			bounds_[d].get_vector( curr_pos_d, best_pos_d );

		Value const vec_to_global_best =
			bounds_[d].get_vector( curr_pos_d, global_best_position[d] );

		Value const unbounded_velocity =
			W * particle.get_current_speed()[ d ] +
			rand_01() * c1 * vec_to_local_best +
			rand_01() * c2 * vec_to_global_best;

		//apply velocity limits
		constexpr Value v_limit_frac = 0.2;
		Value const max_v = v_limit_frac * bounds_[ d ].span();
		Value const min_v = -1.0 * max_v;

		Value const new_speed =
			std::min( std::max( unbounded_velocity, min_v ), max_v );
		particle.set_current_speed( d, new_speed );

		Value const new_pos = bounds_[ d ].apply( curr_pos_d + new_speed );
		particle.set_current_position( d, new_pos );
	}

}


} //namespace

// Reading:
// https://stackoverflow.com/questions/11852577/my-particle-swarm-optimization-code-generates-different-answers-in-c-and-matla
