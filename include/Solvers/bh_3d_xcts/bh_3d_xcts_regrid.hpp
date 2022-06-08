/*
 * Copyright 2022
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 * L. Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
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
#include "kadath_adapted_bh.hpp"
#include "bco_utilities.hpp"
#include "Configurator/config_bco.hpp"
#include "Solvers/co_solver_utils.hpp"
#include <math.h>
#include <sstream>
#include <cmath>
#include <memory>

using namespace Kadath;

template<typename config_t>
int bh_3d_xcts_regrid(config_t& bconfig, const int new_res, std::string outputfile) {
  int exit_status = EXIT_SUCCESS;

  std::string in_spacefile = bconfig.space_filename();
  
  if(!fs::exists(in_spacefile)) {
    std::cerr << "File: " << in_spacefile << " not found.\n\n";
    std::_Exit(EXIT_FAILURE);
  }

  using space_t = Space_adapted_bh;
  FILE *ff1 = fopen(in_spacefile.c_str(), "r");  
  space_t old_space(ff1) ;
  Scalar old_conf (old_space, ff1) ;
	Scalar old_lapse(old_space, ff1) ;
	Vector old_shift(old_space, ff1) ;
  fclose(ff1) ;

	// FIXME not sure if it's only about oddness...
  if(((int)new_res % 2) == 0){
    std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)" << std::endl;
    std::_Exit(EXIT_FAILURE);
  }
  
  // Make sure the domains have adequate space between them based on
  // the previous solution
  bconfig.set(RMID) = bco_utils::get_radius(old_space.get_domain(1), OUTER_BC) ;
  
  // estimate how small the inner radius should be based on relation
  // between conformal factor and numerical radius.
  // see https://arxiv.org/pdf/0805.4192, eq(64)
  double conf_inner = bco_utils::get_boundary_val(2, old_conf, INNER_BC);
  double conf_i_sq  = conf_inner * conf_inner;
  double est_r_div2 = bconfig(MCH) / conf_i_sq;
  bconfig.set(RIN)  =  est_r_div2;
  bconfig.set(ROUT) = 4. * est_r_div2 * 2. ;

  bconfig.set(BCO_RES) = new_res;
  bconfig.set_filename(outputfile);

  // Set this to use radii in setup_co
  bconfig.control(USE_CONFIG_VARS) = true;

  // Generate new space and fields  
  setup_co<BH>(bconfig);
  
  // Unset to prevent possible issues on ID reuse
  bconfig.control(USE_CONFIG_VARS) = false;

  // Read in new fields
  auto new_spacefile = bconfig.space_filename();
  FILE *ff2 = fopen(new_spacefile.c_str(), "r");  
  space_t space(ff2) ;
  Scalar conf (space, ff2) ;
	Scalar lapse(space, ff2) ;
	Vector shift(space, ff2) ;
  fclose(ff2) ;

  std::cout << "Resolution of old space: ";
  bco_utils::print_constant_space_resolution(old_space);

  std::cout << "Resolution of new space: ";
  bco_utils::print_constant_space_resolution(space);

  std::cout << "\nold bounds:" << std::endl;
  bco_utils::print_bounds_from_space(old_space);  
	
  std::cout << "New bounds:" << std::endl;
  bco_utils::print_bounds_from_space(space);
  std::cout << endl;

  // needed in some cases, since the interpolation can go crazy
  const Domain_shell_outer_homothetic* old_outer_homothetic = 
    dynamic_cast<const Domain_shell_outer_homothetic*>(old_space.get_domain(1));
  
  //import fields
  bco_utils::update_adapted_field(old_conf , 2, 1, old_outer_homothetic, OUTER_BC);
  bco_utils::update_adapted_field(old_lapse , 2, 1, old_outer_homothetic, OUTER_BC);
  for(int i = 1; i <= 3; ++i)
    bco_utils::update_adapted_field(old_shift.set(i), 2, 1, old_outer_homothetic, OUTER_BC);

  conf.import(old_conf);
  lapse.import(old_lapse);

  shift.set(1).import(old_shift.set(1));
  shift.set(2).import(old_shift.set(2));
  shift.set(3).import(old_shift.set(3));
  //end import fields

  //reset fields to zero inside excised region
  for(auto& d : {0,1}) {
    lapse.set_domain(d).annule_hard();
    conf.set_domain(d).annule_hard();
    for(int i = 1; i <= 3; ++i) {
      shift.set(i).set_domain(d).annule_hard();
    }
  }
  lapse.std_base();
  conf.std_base();
  shift.std_base();
  
  bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
  return exit_status;
}

