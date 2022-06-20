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
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
#include "coord_fields.hpp"
#include "EOS/EOS.hh"
#include "name_tools.hpp"
#include <cstdlib>
#include <string>
#include <filesystem>
namespace fs = std::filesystem;

inline std::string get_cos_path() {
  const std::string home_kadath{std::getenv("HOME_KADATH")};
  std::string central_abs{home_kadath+"/COs/"};
  return central_abs;
}

#define FORMAT std::setw(13) << std::left << std::showpos 
using namespace Kadath;
const int RELOAD_FILE = 2;
const int RUN_BOOST = 3;

template<typename config_t, typename space_t>
class Solver {
  public:
  using cfgen_t  = CoordFields<space_t>;
  using cfary_t  = std::array<std::optional<Vector>, NUM_VECTORS>;
  using base_config_t = std::decay_t<config_t>;
  using base_space_t = std::decay_t<space_t>;
 
  protected:
  space_t& space;
	config_t& bconfig;
  Base_tensor& basis;
  cfgen_t cfields;
  cfary_t coord_vectors;
  const int ndom;

  public:
  Solver(config_t& config_in, space_t& space_in, Base_tensor& base_in) 
    : space(space_in), basis(base_in), bconfig(config_in), 
        cfields(space), ndom(space_in.get_nbr_domains()) {}
  virtual void print_diagnostics(const System_of_eqs& syst, const int ite, 
    const double conv) const = 0;
  virtual std::string converged_filename(const std::string& stage) const = 0;
  virtual ~Solver() = default;

  void check_max_iter_exceeded(const int& rank, const int& ite, const double& conv) const {
    bool exceeded = false;
    if(ite > bconfig.seq_setting(MAX_ITER) && \
       conv < 10. * bconfig.seq_setting(PREC)) {
      if(rank == 0)
        std::cout << "Max iterations exceeded at precison, " << conv
                  << "\nFinishing since precision < 10. * PREC....\n"
                  << "Running at higher resolution may help.\n";
      std::_Exit(EXIT_FAILURE);
    }
    else if(ite > bconfig.seq_setting(MAX_ITER)) {
      if(rank == 0)
        std::cout << "Max iterations exceeded at precison, " << conv;
      std::_Exit(EXIT_FAILURE);
    }
  }

  bool solution_exists(std::string last_stage="") {
    bool exists = false;
    std::string prev_name{converged_filename(last_stage)};
    std::string prev_abs{bconfig.config_outputdir()+prev_name};
    
    const std::string home_kadath{std::getenv("HOME_KADATH")};
    std::string central_abs{home_kadath+"/COs/"+prev_name};
    
    auto check_and_update_config = [&](auto p) {
      exists = true;
      base_config_t old_solution(p+".info");
      if(bconfig.set(NSHELLS) != old_solution.set(NSHELLS))
        return false;
      
      // make sure we copy stages, controls, and settings over
      for(size_t idx = 0; idx < NUM_STAGES; ++idx)
        old_solution.set_stage(idx) = bconfig.set_stage(idx);
      for(size_t idx = 0; idx < NUM_CONTROLS; ++idx)
        old_solution.control(idx) = bconfig.control(idx);
      for(size_t idx = 0; idx < NUM_SEQ_SETTINGS; ++idx)
        old_solution.seq_setting(idx) = bconfig.seq_setting(idx);
      
      bconfig = old_solution;
      return exists;
    };

    if(fs::exists(prev_abs+".info") && fs::exists(prev_abs+".dat")){
      exists = check_and_update_config(prev_abs);
    }
    else if(bconfig.control(SAVE_COS) && fs::exists(central_abs+".info") && fs::exists(central_abs+".dat")) {
      exists = check_and_update_config(central_abs);
    }
    
    return exists;
  }

  template<typename... bco_idx>
  std::string extract_eos_name(bco_idx... bco) const {
    const std::string eos_file_abs = bconfig.template eos<std::string>(EOSFILE, bco...);
    const std::string eos_file = extract_filename(eos_file_abs);
    return eos_file.substr(0, eos_file.find("."));
  }
};

