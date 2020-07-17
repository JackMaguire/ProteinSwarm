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

  void set_lower_bound( int const dim, double const lb ){
    bounds_[ dim ].lower_bound = lb;
  }

  void set_all_upper_bounds( double const ub ){
    for( Bounds & b : bounds_ ){
      b.upper_bound = ub;
    }
  }

  void set_upper_bound( int const dim, double const ub ){
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

    bounds_[ dim ].type = t;
  }

  ProteinSwarm make_pso( int const n_particles ){
    protein_swarm::ProteinSwarm swarm( n_particles, bounds_ );
    return swarm;
  }

  std::vector< Bounds > bounds_;
};

//https://stackoverflow.com/questions/10701514/how-to-return-numpy-array-from-boostpython
template< typename T >
np::ndarray
vec_to_np( std::vector< T > const & v ){
  Py_intptr_t shape[1] = { (Py_intptr_t) v.size() };
  np::ndarray result = np::zeros(1, shape, np::dtype::get_builtin< T >());
  std::copy( v.begin(), v.end(), reinterpret_cast<T*>(result.get_data() ) );
  return result;
}

np::ndarray
extract_data( Sample sample ) {
  return vec_to_np( sample.value );
}

np::ndarray
get_best_position( ProteinSwarm const & optimizer ) {
  //I don't think pass-by-ref works with python

  return vec_to_np( optimizer.get_global_best_position() );
}

BOOST_PYTHON_MODULE( proteinswarm )
{
  using namespace boost::python;
  Py_Initialize();
  np::initialize();

  void (ProteinSwarm::*tell)(Sample const &,Value) = &ProteinSwarm::tell;
  void (ProteinSwarm::*tell_sampleinfo)(SampleInfo, Value) = &ProteinSwarm::tell;
  void (ProteinSwarm::*tell_particle_id)(uint, Value) = &ProteinSwarm::tell;

  def( "extract_data", &extract_data );
  def( "get_best_position", &get_best_position );
  
  class_< Sample >( "Sample" )
    .def( "get_particle_id", &Sample::get_particle_id );

  
  class_< ProteinSwarm >( "ProteinSwarm", init<uint,std::vector< Bounds > const &,InitialSamplingMethod>() )
    .def( "reset", &ProteinSwarm::reset )
    .def( "ask", &ProteinSwarm::ask )
    .def( "tell", tell )
    .def( "tell_sampleinfo", tell_sampleinfo )
    .def( "get_global_best_score", &ProteinSwarm::get_global_best_score )
    .def( "set_span_coeff", &ProteinSwarm::set_span_coeff )
    .def( "set_k1", &ProteinSwarm::set_k1 )
    .def( "set_k2", &ProteinSwarm::set_k2 )
    .def( "set_c1", &ProteinSwarm::set_c1 )
    .def( "set_c2", &ProteinSwarm::set_c2 )
    .def( "set_v_limit_m", &ProteinSwarm::set_v_limit_m )
    .def( "set_v_limit_b", &ProteinSwarm::set_v_limit_b )
    .def( "parameters_are_reasonable", &ProteinSwarm::parameters_are_reasonable )
    .def( "get_n_particles", &ProteinSwarm::get_n_particles );

  class_< ProteinSwarmBounds >( "ProteinSwarmBounds", init<int>() )
    .def( "set_all_lower_bounds", &ProteinSwarmBounds::set_all_lower_bounds )
    .def( "set_lower_bound", &ProteinSwarmBounds::set_lower_bound )
    .def( "set_all_upper_bounds", &ProteinSwarmBounds::set_all_upper_bounds )
    .def( "set_upper_bound", &ProteinSwarmBounds::set_upper_bound )
    .def( "set_all_bounds_type", &ProteinSwarmBounds::set_all_bounds_type )
    .def( "set_bounds_type", &ProteinSwarmBounds::set_bounds_type )
    .def( "make_pso", &ProteinSwarmBounds::make_pso );

}
