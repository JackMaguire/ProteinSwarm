#include "../include/protein_swarm.hh"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <vector>
#include <cctype> //tolower
#include <stdexcept>

using namespace protein_swarm;

using namespace boost;
using namespace boost::python;
using namespace boost::python::numpy;

namespace p = boost::python;
namespace np = boost::python::numpy;

struct ProteinSwarmBounds {
  ProteinSwarmBounds( int ndim ) {
    bounds_.resize( ndim );
  }

  void set_all_lower_bounds( double const lb ){
    for( Bounds & b : bounds_ ){
      b.lower_bound = lb;
    }
  }

  void set_all_lower_bounds( int const dim, double const lb ){
    bounds_[ dim ].lower_bound = lb;
  }

  void set_all_upper_bounds( double const ub ){
    for( Bounds & b : bounds_ ){
      b.upper_bound = ub;
    }
  }

  void set_all_upper_bounds( int const dim, double const ub ){
    bounds_[ dim ].upper_bound = ub;
  }

  void set_all_bounds_type(
    std::string type //pass-by-value on purpose
  ){
    BoundsType t;
    if( type == "standard" ) t = BoundsType::STANDARD;
    else if( type == "pacman" ) t = BoundsType::PACMAN;
    else throw std::runtime_error( type + " is not an option. Must be 'standard' or 'pacman'" );

    for( Bounds & b : bounds_ ){
      b.type = t;
    }
  }

  void set_bounds_type(
    int dim,
    std::string type //pass-by-value on purpose
  ){
    BoundsType t;
    if( type == "standard" ) t = BoundsType::STANDARD;
    else if( type == "pacman" ) t = BoundsType::PACMAN;
    else throw std::runtime_error( type + " is not an option. Must be 'standard' or 'pacman'" );

    bounds_[ dim ] = t;
  }

  ProteinSwarm make_pso( int const n_particles ){
    protein_swarm::ProteinSwarm swarm( n_particles, bounds_.size(), bounds_ );
    return swarm;
  }

  std::vector< Bounds > bounds_;
};

BOOST_PYTHON_MODULE( proteinswarm )
{
  using namespace boost::python;
  Py_Initialize();
  np::initialize();

  void (ProteinSwarm::*tell_sample)(Sample const &,Value) = &ProteinSwarm::tell;
  void (ProteinSwarm::*tell_sampleinfo)(SampleInfo, Value) = &ProteinSwarm::tell;

  class_< ProteinSwarm >( "ProteinSwarm" )
    .def( "reset", &ProteinSwarmBounds::reset )
    .def( "ask", &ProteinSwarmBounds::ask )
    .def( "tell_sample", tell_sample )
    .def( "tell_sampleinfo", tell_sampleinfo )
    .def( "get_n_particles", &ProteinSwarmBounds::get_n_particles );

  class_< ProteinSwarmBounds >( "ProteinSwarmBounds" )
    .def( "set_all_lower_bounds", &ProteinSwarmBounds::set_all_lower_bounds )
    .def( "set_lower_bounds", &ProteinSwarmBounds::set_lower_bounds )
    .def( "set_all_upper_bounds", &ProteinSwarmBounds::set_all_upper_bounds )
    .def( "set_upper_bounds", &ProteinSwarmBounds::set_upper_bounds )
    .def( "set_all_bounds_type", &ProteinSwarmBounds::set_all_bounds_type )
    .def( "set_bounds_type", &ProteinSwarmBounds::set_bounds_type )
    .def( "make_pso", &ProteinSwarmBounds::make_pso );

}
