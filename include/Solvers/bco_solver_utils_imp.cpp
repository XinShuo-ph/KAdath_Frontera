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
#include "solvers.hpp"
#include "mpi.h"
#include "bco_utilities.hpp"
#include "bco_solver_utils.hpp"
#include "bns_xcts/bns_xcts_solver.hpp"
#include "ns_3d_xcts/ns_3d_xcts_solver.hpp"
#include "bh_3d_xcts/bh_3d_xcts_solver.hpp"
#include <iostream>

template<typename config_t>
std::string solve_NS_from_binary(config_t& bconfig, const size_t bco) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::string output_path = (bconfig.control(SAVE_COS)) ? get_cos_path() : "./COs";
  fs::create_directory(output_path);
  
  kadath_config_boost<BCO_NS_INFO> nsconfig;
  nsconfig.set_defaults();

  // Tells the NS driver to initialize the numerical space and fields
  nsconfig.control(SEQUENCES) = true;
  
  // copy parameters from binary configuration
  for(int i = 0; i < NUM_BCO_PARAMS; ++i) 
    nsconfig.set(i) = bconfig.set(i, bco);

  for(int i = 0; i < NUM_EOS_PARAMS; ++i) 
    nsconfig.set_eos(i) = bconfig.set_eos(i, bco);

  for(int i = 0; i < NUM_CONTROLS; ++i)
    nsconfig.control(i) = bconfig.control(i);

  for(int i = 0; i < NUM_SEQ_SETTINGS; ++i)
    nsconfig.seq_setting(i) = bconfig.seq_setting(i);
  nsconfig.seq_setting(INIT_RES) = 9;
  nsconfig.set_outputdir(output_path);
  
  if(bconfig.control(USE_BOOSTED_CO))
    nsconfig.set_stage(BIN_BOOST) = true;
  
  const int ninshells = bconfig(NINSHELLS, bco);
  const int nshells = bconfig(NSHELLS, bco);

  nsconfig(NINSHELLS) = 0;
  nsconfig(NSHELLS) = (bconfig.control(CO_USE_SHELLS)) ?
    nshells : 0;
  
  int err = ns_3d_xcts_driver(nsconfig, output_path, bconfig);
  MPI_Barrier(MPI_COMM_WORLD);

  // update binary parameters
  for(int i = 0; i < NUM_BCO_PARAMS; ++i) 
    bconfig.set(i, bco) = nsconfig.set(i) ;
  // put shells back for the binary
  bconfig.set(NINSHELLS, bco) = ninshells;
  bconfig.set(NSHELLS, bco) = nshells;
  
  return nsconfig.config_outputdir()+nsconfig.config_filename();
}

template<typename config_t>
std::string solve_BH_from_binary(config_t& bconfig, const size_t bco) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::string output_path = (bconfig.control(SAVE_COS)) ? get_cos_path() : "./COs";
  fs::create_directory(output_path);
  
  kadath_config_boost<BCO_BH_INFO> bhconfig;
  bhconfig.set_defaults();
  
  // Tells the NS driver to initialize the numerical space and fields
  bhconfig.control(SEQUENCES) = true;
  
  // copy parameters from binary configuration
  for(int i = 0; i < NUM_BCO_PARAMS; ++i) 
    bhconfig.set(i) = bconfig.set(i, bco);
  
  for(int i = 0; i < NUM_CONTROLS; ++i)
    bhconfig.control(i) = bconfig.control(i);

  for(int i = 0; i < NUM_SEQ_SETTINGS; ++i)
    bhconfig.seq_setting(i) = bconfig.seq_setting(i);
  bhconfig.seq_setting(INIT_RES) = 9;
  bhconfig.set_filename("initbh");
  bhconfig.set_outputdir(output_path);
  if(bconfig.control(USE_BOOSTED_CO))
    bhconfig.set_stage(BIN_BOOST) = true;
  
  const int nshells = bconfig(NSHELLS, bco);
  bhconfig(NSHELLS) = (bconfig.control(CO_USE_SHELLS)) ?
    nshells : 0;
  
  int err = bh_3d_xcts_driver(bhconfig, output_path, bconfig, bco);
  MPI_Barrier(MPI_COMM_WORLD);

  // update binary parameters
  for(int i = 0; i < NUM_BCO_PARAMS; ++i) 
    bconfig.set(i, bco) = bhconfig.set(i) ;
  bconfig.set(NSHELLS, bco) = nshells;
 
  MPI_Barrier(MPI_COMM_WORLD);
  return bhconfig.config_outputdir()+bhconfig.config_filename();
}

inline
void check_dist(double dist, double M1, double M2, double garbage_factor) {
  auto M = M1 + M2;
  auto const garbage_dist = garbage_factor * M;
  auto const recommended_dist = 8. * M;
  if(dist <= garbage_dist) {
    std::cerr << "Distance is set to (" << dist
              << ") which will not give results. \nSet to ("
              << recommended_dist << ") for something reasonable.\n";
    std::_Exit(EXIT_FAILURE);
  }
}
