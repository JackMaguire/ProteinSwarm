// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-

#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <assert.h>
#include <limits>

//#include <array>
#include <vector>
#include <queue>

/*
  CURRENT ASSUMPTIONS
  - N Partcles >= N CPUs. There will never be more than one CPU working on the same particle at the same time
*/

#define rand_01 ((float)rand() / (float)RAND_MAX)

namespace protein_swarm {

using uint = unsigned int;
using ulong = unsigned long int;

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

template< uint NDIM, typename Value = double >
class ParticleInitialzer {
public:
  virtual ~ParticleInitialzer(){};

  //This will be called many times for a single run, once for each particle
  virtual
  std::vector< Value >
  initialize_one_particle(
    uint n_total_particles, //How many particles are in the system?
    uint particle_id, //which particle are we considering now? Ranging from [ 0, n_total_particles-1 ] inclusive
    std::vector< Value > const & lower_bounds,//these are the bounds, which the user probably already knows anyways
    std::vector< Value > const & upper_bounds
  ) const = 0; //returns the starting values for each dimension for this particle

};

struct SampleInfo {
  uint particle;
};

//This is one case where array of struct appears to be better than struct of array
template< typename Value = double >
class ParticleInfo {
private:
	Value current_score_ = std::numeric_limits< Value >::max();
	std::vector< Value > current_position_; //Keep empty until it is set
	std::vector< Value > current_speed_;

	Value best_score_ = std::numeric_limits< Value >::max();
	std::vector< Value > best_position_; //Keep empty until it is set

	bool current_score_is_being_evaluated_ = true;
	//This is true when the current_position_ has been set to the new value but the current_score_ has not been calculated yet. current_score_is_being_evaluated_ == true means that the job is currently in flight

public:
	//ParticleInfo() = default;
	ParticleInfo( ParticleInfo const & ) = default;
	ParticleInfo( ParticleInfo && ) = default;

	ParticleInfo(){
		//current_position_.resize( ndim ); //keep positions empty until they are set
		//current_speed_.resize( ndim ); //keep empty until it is set
		//best_position_.resize( ndim ); //keep positions empty until they are set
	}

	~ParticleInfo() = default;

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

	void
	update_to_new_position(
		std::vector< Value > const & global_best_position,
		std::vector< Value > const & lower_bounds,
		std::vector< Value > const & upper_bounds
	){
		assert( ! current_position_.empty() ); //Did you call initialize()?
		assert( ! current_speed_.empty() ); //Did you call initialize()?

		//sample a new position
		uint const ndim = current_position_.size();
		for( uint d = 0; d < ndim; ++d ){
			Value velocity = 
		}
	}

public: //accessors
	Value get_current_score() const {
		return current_score_;
	}

	void set_current_score( Value v ){
		current_score_ = v;
		current_score_is_being_evaluated_ = false;
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

	std::vector< Value > const & get_current_speed() const {
		return current_speed_;
	}

	void set_current_speed( std::vector< Value > const & v ){
		current_speed_ = v;
	}

	void set_current_speed( std::vector< Value > && v ){
		current_speed_ = v;
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

};

//THIS IS NOT THREADSAFE
template< typename Value = double >
class ProteinSwarm {

  struct Sample {
    SampleInfo info;
    std::vector< Value > value;
  };

public:
  ProteinSwarm(
    uint const n_particles,
    std::vector< Value > const & lower_bounds,
    std::vector< Value > const & upper_bounds,
    InitialSamplingMethod const sampling_method = UNIFORM
  ):
    n_particles_( n_particles ),
    lower_bounds_( lower_bounds ),
    upper_bounds_( upper_bounds )
  {
    switch( sampling_method ){
      //TODO
    }
    initialize();
  }

  void reset(){
    //TODO
  }

	Sample ask(){
		assert( ! particle_queue_.empty() );
		uint const particle = particle_queue_.front();
		particle_queue_.pop();
		//TODO
		++n_asks_;
	}

  void tell(
    Sample const & sample,
    Value const score
  ){
    tell( sample.info, score );
  }

  void tell(
    SampleInfo const info,
    Value const score
  ){
    //TODO
		++n_tells;
  }

public://getters and setters
  uint get_n_particles() const {
    return n_particles_;
  }

  void set_n_particles( uint setting ){
    n_particles_ = setting;
  }

protected:
  void initialize(
    ParticleInitialzer const & initializer
  ){
    particles_.resize( n_particles_ );
    for( uint i = 0; i < n_particles_; ++i ){
			particle_queue_.push( i );

			particles_[ i ].current_position

      //current_positions_[ i ] = initializer.initialize_one_particle(	n_particles_, i, lower_bounds_, upper_bounds_ );
    }
    //TODO
  }

private:
  uint n_particles_;
	uint ndim_;

	//ParticleInfo holds location, history, and score info for each particle
	std::vector< ParticleInfo > particles_;
  uint index_of_global_best_ = 0;

	//The particle queue holds a list of particles that are ready to be sampled
	std::queue< uint > particle_queue_;

  //bool bounded_ = false;
  std::vector< Value > lower_bounds_;
  std::vector< Value > upper_bounds_;

  //uint njobs_submitted_ = 0;
  uint n_asks_ = 0;
	uint n_tells_ = 0;

};

}

const int numofdims = 30;
const int numofparticles = 50;

using namespace std;

void fitnessfunc(float X[numofparticles][numofdims], float fitnesses[numofparticles])
{
  memset(fitnesses, 0, sizeof (float) * numofparticles);
  for(int i = 0; i < numofparticles; i++)
    {
      for(int j = 0; j < numofdims; j++)
        {
	  fitnesses[i] += (pow(X[i][j], 2));
        }
    }
}

float mean(float inputval[], int vallength)
{
  int addvalue = 0;
  for(int i = 0; i < vallength; i++)
    {
      addvalue += inputval[i];
    }
  return (float)(addvalue / vallength);
}

void PSO(int numofiterations, float c1, float c2,
  float Xmin[numofdims], float Xmax[numofdims], float initialpop[numofparticles][numofdims],
  float worsts[], float meanfits[], float bests[], float *gbestfit, float gbest[numofdims])
{
  float V[numofparticles][numofdims] = {0};
  float X[numofparticles][numofdims];
  float Vmax[numofdims];
  float Vmin[numofdims];
  float pbests[numofparticles][numofdims];
  float pbestfits[numofparticles];
  float fitnesses[numofparticles];
  float w;
  float minfit;
  int   minfitidx;

  memcpy(X, initialpop, sizeof(float) * numofparticles * numofdims);
  fitnessfunc(X, fitnesses);
  minfit = *min_element(fitnesses, fitnesses + numofparticles);
  minfitidx = min_element(fitnesses, fitnesses + numofparticles) - fitnesses;
  *gbestfit = minfit;
  memcpy(gbest, X[minfitidx], sizeof(float) * numofdims);

  for(int i = 0; i < numofdims; i++)
    {
      Vmax[i] = 0.2 * (Xmax[i] - Xmin[i]);
      Vmin[i] = -Vmax[i];
    }

  for(int t = 0; t < 1000; t++)
    {
      w = 0.9 - ( (0.7 * t) / numofiterations);

      for(int i = 0; i < numofparticles; i++)
        {
	  if(fitnesses[i] < pbestfits[i])
            {
	      pbestfits[i] = fitnesses[i];
	      memcpy(pbests[i], X[i], sizeof(float) * numofdims);
            }
        }
      for(int i = 0; i < numofparticles; i++)
        {
	  for(int j = 0; j < numofdims; j++)
            {
	      V[i][j] = min(max((w * V[i][j] + rand_01 * c1 * (pbests[i][j] - X[i][j])
		    + rand_01 * c2 * (gbest[j] - X[i][j])), Vmin[j]), Vmax[j]);
	      X[i][j] = min(max((X[i][j] + V[i][j]), Xmin[j]), Xmax[j]);
            }
        }

      fitnessfunc(X, fitnesses);
      minfit = *min_element(fitnesses, fitnesses + numofparticles);
      minfitidx = min_element(fitnesses, fitnesses + numofparticles) - fitnesses;
      if(minfit < *gbestfit)
        {
	  *gbestfit = minfit;
	  memcpy(gbest, X[minfitidx], sizeof(float) * numofdims);
        }

      worsts[t] = *max_element(fitnesses, fitnesses + numofparticles);
      bests[t] = *gbestfit;
      meanfits[t] = mean(fitnesses, numofparticles);
    }


}

int main()
{
  time_t t;
  srand((unsigned) time(&t));

  float xmin[30], xmax[30];
  float initpop[50][30];
  float worsts[1000], bests[1000];
  float meanfits[1000];
  float gbestfit;
  float gbest[30];
  for(int i = 0; i < 30; i++)
    {
      xmax[i] = 100;
      xmin[i] = -100;
    }
  for(int i = 0; i < 50; i++)
    for(int j = 0; j < 30; j++)
      {
	initpop[i][j] = rand() % (100 + 100 + 1) - 100;
      }

  PSO(1000, 2, 2, xmin, xmax, initpop, worsts, meanfits, bests, &gbestfit, gbest);

  cout<<"fitness: "<<gbestfit<<endl;
  return 0;
}

// Reading:
// https://stackoverflow.com/questions/11852577/my-particle-swarm-optimization-code-generates-different-answers-in-c-and-matla
