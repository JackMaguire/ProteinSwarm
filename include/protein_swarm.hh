#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>

#include <array>
#include <vector>

#define rand_01 ((float)rand() / (float)RAND_MAX)

namespace protein_swarm {

using uint = unsigned int;
using ulong = unsigned long int;

enum class InitialSamplingMethod {

  //every particle randomly chooses its starting position from [ lower_bound, upper_bound ] uniformly
  UNIFORM,

  //every particle randomly chooses its starting position from [ -1, 1 ] uniformly
  UNIFORM_NORM,

  //every particle randomly chooses its starting position from [ 0, 1 ] uniformly
  UNIFORM_NORM_0,

  //every particle randomly chooses its starting position by sampling from a gaussian distribution centered at lower_bound + (upper_bound-lower_bound)/2
  //with a std deviation of (upper_bound-lower_bound)/6 and clipped at [ lower_bound, upper_bound ]
  GAUSSIAN,

  //every particle randomly chooses its starting position by sampling from a gaussian distribution centered at 0
  //with a std deviation of (upper_bound-lower_bound)/6 and clipped at [ lower_bound, upper_bound ]
  GAUSSIAN_0,

  //every particle randomly chooses its starting position by sampling from a gaussian distribution centered at 0
  //with a std deviation of 1 and clipped at [ lower_bound, upper_bound ]
  GAUSSIAN_NORM,

  //particles attempt to be evenly spaced out
  GRID,

  //use user-provided intialization by calling set_particle_initializer
  CUSTOM
};

template< uint NDIM, typename Value = double >
class ParticleInitialzer {
public:
  virtual ~ParticleInitialzer(){};

  //This will be called many times for a single run, once for each particle
  virtual
  std::array< Value, NDIM >
  initialize_one_particle(
    uint n_total_particles, //How many particles are in the system?
    uint particle_id, //which particle are we considering now? Ranging from [ 0, n_total_particles-1 ] inclusive
    std::array< Value, NDIM > const & lower_bounds,//these are the bounds, which the user probably already knows anyways
    std::array< Value, NDIM > const & upper_bounds
  ) const = 0; //returns the starting values for each dimension for this particle

};

struct SampleInfo {
  uint particle;
};


//THIS IS NOT THREADSAFE
template< uint NDIM, typename Value = double >
class ProteinSwarm {

struct Sample {
  SampleInfo info;
  std::array< Value, NDIM > value;
};

public:
  ProteinSwarm(
    uint const n_particles,
    std::array< Value, NDIM > const & lower_bounds,
    std::array< Value, NDIM > const & upper_bounds,
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
    uint const particle = n_asks_general_ % n_particles_;
    ++n_asks_general_;
    return ask( particle );
  }

  Sample ask(
    uint const particle
  ){
    //TODO

    ++n_asks_total_;
    ++n_asks_in_flight_;
  }

  void tell(
    Sample const & sample,
    double const score
  ){
    tell( sample.info, score );
  }

  void tell(
    SampleInfo const info,
    double const score
  ){
    //TODO
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
    current_positions_.resize( n_particles_ );
    for( uint i = 0; i < n_particles_; ++i ){
      //auto == std::array< Value, NDIM >
      current_positions_[ i ] = initializer.initialize_one_particle(
	n_particles_, i, lower_bounds_, upper_bounds_ );
    }
    //TODO
  }

private:
  //bool bounded_ = false;
  std::array< Value, NDIM > lower_bounds_;
  std::array< Value, NDIM > upper_bounds_;

  //uint njobs_submitted_ = 0;
  uint n_asks_general_ = 0;
  uint n_asks_total_ = 0;
  uint n_asks_in_flight_ = 0;

  uint n_particles_;
  std::vector< std::array< Value, NDIM > > current_positions_;
  std::vector< double > scores_for_current_positions_;

  std::vector< std::array< Value, NDIM > > best_positions_;
  std::vector< double > scores_for_best_positions_;
  uint index_of_global_best_ = 0;
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
