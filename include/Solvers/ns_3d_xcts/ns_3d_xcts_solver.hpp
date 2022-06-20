#pragma once
#include "Solvers/solvers.hpp"
using namespace Kadath;

template<class eos_t, typename config_t, typename space_t = Kadath::Space_spheric_adapted>
class ns_3d_xcts_solver : Solver<config_t, space_t> {
  public:
  using typename Solver<config_t, space_t>::base_config_t;
  using typename Solver<config_t, space_t>::base_space_t;

  private:
  Scalar& conf;
  Scalar& lapse;
  Scalar& logh;
  Vector& shift;
  Metric_flat fmet;

  // Specify base class members used to avoid this->
  using Solver<config_t, space_t>::space;
  using Solver<config_t, space_t>::bconfig;
  using Solver<config_t, space_t>::basis;
  using Solver<config_t, space_t>::cfields;
  using Solver<config_t, space_t>::coord_vectors;
  using Solver<config_t, space_t>::ndom;
  using Solver<config_t, space_t>::check_max_iter_exceeded;
  using Solver<config_t, space_t>::solution_exists;
  using Solver<config_t, space_t>::extract_eos_name;

  public:
  // solver is not trivially constructable since Kadath containers are not
  // trivially constructable
  ns_3d_xcts_solver() = delete;

  ns_3d_xcts_solver(config_t& config_in, space_t& space_in, Base_tensor& base_in,  
    Scalar& conf_in, Scalar& lapse_in, Scalar& logh_in, Vector& shift_in);
  
  // syst always requires the same initialization for the stages
  void syst_init(System_of_eqs& syst);
  
  // diagnostics at runtime
  void print_diagnostics_norot(const System_of_eqs& syst, 
    const int ite, const double conv) const;
  void print_diagnostics(const System_of_eqs& syst, 
    const int  ite = 0, const double conv = 0) const override;
  
  std::string converged_filename(const std::string& stage="") const override;
  //void update_stages(config_t& old_config);
  
  // solve driver
  int solve();

  // solver stages
  int norot_stage(bool fixed = false);
  int uniform_rot_stage();
  
  /**
   * binary_boost_stage
   *
   * based on an input binary Configurator file, we boost the BH accordingly
   *
   * @param[input] binconfig: binary Configurator file
   * @param[input] bco: index of BCO - needed to determine coordinate shift based on BCO location in binary space
   */
  int binary_boost_stage(kadath_config_boost<BIN_INFO>& binconfig, const size_t bco);

  // Update bconfig(HC) and bconfig(NC)
  void update_config_quantities(const double& loghc);
 
//  template<class eos_t>
  //bool useful_checkpoint_exists();
};

/**
 * ns_3d_xcts_driver
 *
 * Control computation of 3D NS from setup config/dat file combination
 * filename is pulled from bconfig which should pair with <filename>.dat
 *
 * This code exists to allow the required  NS' to be computed during a BNS solver
 * while also minimizing the code in an isolated solver.
 *
 * Also, KADATH at the time of writing, was strict on not allowing trivial
 * construction of base types (Base_tensor, Scalar, Tensor, Space, etc) which
 * means a driver is required to populate the related Solver class.
 *
 * @tparam[int] bconfig - Configurator file
 * @return Success/failure
 */
template<typename config_t>
inline int ns_3d_xcts_driver (config_t& bconfig, std::string outputdir, 
 kadath_config_boost<BIN_INFO> binconfig = kadath_config_boost<BIN_INFO>{}, 
 const size_t bco = BCO1);


#include "ns_3d_xcts_solver_imp.cpp"
#include "ns_3d_xcts_stages.cpp"
