/*
 * Copyright 2021
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
 *
 */
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "kadath.hpp"
#include "kadath_adapted_polar.hpp"
#include "kadath_adapted.hpp"
#include <fstream>
#include <string>
#include <type_traits>
#include "EOS/EOS.hh"

using namespace Kadath;
using namespace Kadath::Margherita;

template<typename Config>
int norot_3dsetup (Config bconfig);

int main(int argc, char** argv) {
  // class to load the configuration
  kadath_config_boost<BCO_NS_INFO> bconfig;

  // check if configuration is given
  // or if a standard setup should be generated
  if (argc < 2) {
    std::cout << "Using default Togashi setup.\n";
    bconfig.set_defaults();
  }
  else {
    std::string input_filename = std::string{argv[1]};
    std::cout << "Creating setup from: " << input_filename << "\n";
    
    bconfig.set_filename(input_filename);
    bconfig.open_config();
  }
  
  // set new base filename for new setup
  std::string base_filename = "./initns.info";
  bconfig.set_filename(base_filename);

  /* Loading of the configuration parameters.
   *
   * See $HOME_KADATH/include/Configurator/config_enums.hpp for
   * possible BCO_PARAMS and EOS_PARAMS index names 
   */
  const double h_cut = bconfig.eos<double>(HC);
  const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

  if(eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
    bconfig.set(HC) = EOS<eos_t,PRESSURE>::h_cold__rho(bconfig(NC));
  } else if(eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
                            2000 : bconfig.eos<int>(INTERP_PTS);

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
    bconfig.set(HC) = EOS<eos_t,PRESSURE>::h_cold__rho(bconfig(NC));
  }
  else { 
    std::cerr << eos_type << " is not recognized.\n";
    std::_Exit(EXIT_FAILURE);
  }

  // call setup routine, setting up the actual numerical domains
  int res = norot_3dsetup(bconfig);
  std::cout << bconfig << std::endl;

  return 0;
}

// this function takes a configuration and creates the corresponding numerical setup
template<typename Config>
int norot_3dsetup(Config bconfig) {
  // number of dimensions
  const int dim = bconfig(DIM);

  // collocation point type and number of collocation points per domain
  int type_coloc = CHEB_TYPE;
  Dim_array res(dim);
  res.set(0) = bconfig(BCO_RES);
  res.set(1) = bconfig(BCO_RES);
  res.set(2) = bconfig(BCO_RES)-1;

  // physical center of the system, arbitrary
  Point center(dim);
  for (int i = 1; i <= dim; i++)
    center.set(i) = 0;
  
  if(!bconfig.control(USE_CONFIG_VARS)) {
    bconfig.set(RIN) = 0.5 * bconfig(RMID);
    bconfig.set(ROUT) = 1.5 * bconfig(RMID);
  }

  // domain boundaries, i.e. radii for each domain of the single star space
  int ndom      = 4;
  Array<double> bounds(ndom - 1);
  bounds.set(0) = bconfig(RIN);
  bounds.set(1) = bconfig(RMID);
  bounds.set(2) = bconfig(ROUT);

  // generate a full single star space including compactification to infinity
  Space_spheric_adapted space(type_coloc, center, res, bounds);
  // use a basis of cartesian type
  Base_tensor basis(space, CARTESIAN_BASIS);

  // start to create the fields

  // H = log(h), the logarithm of the enthalpy
  Scalar logh(space);
  // initialize to the central value, which should be picked carefully
  logh.set_domain(0) = std::log(bconfig(HC));
  logh.set_domain(1) = std::log(bconfig(HC));
  // set the coefficients to zero outside of the star
  for (int d = 2; d < ndom; ++d)
    logh.set_domain(d).annule_hard();
  logh.std_base();

  // conformal factor is initialized to one, representing flat space as an initial guess
  Scalar conf(space);
  conf = 1.;
  conf.std_base();

  // same for the lapse
  Scalar lapse(conf);
  lapse.std_base();

  // shift is zero in the beginning and in general for TOV solutions
  Vector shift(space, CON, basis);
  shift.annule_hard();
  shift.std_base();
  // finished create the fields

  // write space and config to files on one processor
  bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh);
  return EXIT_SUCCESS;
}
