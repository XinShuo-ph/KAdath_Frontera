/*
 * Copyright 2022
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once
#include "Solvers/solvers.hpp"
#include "kadath_bin_ns.hpp"
using namespace Kadath;

template<class eos_t, typename config_t, typename space_t = Kadath::Space_bin_ns>
class bns_xcts_solver : Solver<config_t, space_t> {
  public:
  using typename Solver<config_t, space_t>::base_config_t;
  using typename Solver<config_t, space_t>::base_space_t;

  private:
  Scalar& conf;
  Scalar& lapse;
  Scalar& logh;
  Scalar& phi;
  Vector& shift;
  Metric_flat fmet;

  const double xc1;
  const double xc2;
  const double xo{0.};

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
  using Solver<config_t, space_t>::checkpoint;

  public:
  // solver is not trivially constructable since Kadath containers are not
  // trivially constructable
  bns_xcts_solver() = delete;

  bns_xcts_solver(config_t& config_in, space_t& space_in, Base_tensor& base_in,  
    Scalar& conf_in, Scalar& lapse_in, Vector& shift_in, Scalar& logh_in, Scalar& phi_in);
  
  // syst always requires the same initialization for the stages
  void syst_init(System_of_eqs& syst);
  
  // diagnostics at runtime
  void print_diagnostics(const System_of_eqs& syst, 
    const int  ite = 0, const double conv = 0) const override;
  
  std::string converged_filename(const std::string stage="") const override;
  
  void save_to_file() const override {
    bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
  }
  
  // solve driver
  int solve();

  /**
   * hydrostatic_equilibrium stage
   *
   * The binary is solved using the force balance equations
   * in addition to solving the relativisitc Euler equation
   * to obtain a binary in hydrostatic equilibrium.
   */
  int hydrostatic_equilibrium_stage();
  /**
   * hydro_rescaling_stages
   *
   * The binary is solved using fixed orbital velocity
   * therefore the matter scalar fields are simply rescaled
   * based on the fixed baryonic mass of each star.
   *
   * @param[input] stage: Some changes are made based on TOTAL_BC or ECC_RED stage.
   */
  int hydro_rescaling_stages(const size_t stage, std::string stage_text);

  // Update bconfig(HC) and bconfig(NC)
  void update_config_quantities(const double& loghc);
};

template<typename eos_t>
inline void bns_setup_boosted_3d(
  kadath_config_boost<BCO_NS_INFO>& NS1config, kadath_config_boost<BCO_NS_INFO>& NS2config,
  kadath_config_boost<BIN_INFO>& bconfig);

#include "bns_xcts_solver_imp.cpp"
#include "bns_xcts_stages.cpp"
