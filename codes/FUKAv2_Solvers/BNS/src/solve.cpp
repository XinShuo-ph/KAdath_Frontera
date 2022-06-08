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
#include "mpi.h"
#include "Solvers/bns_xcts/bns_xcts_solver.hpp"
#include "Solvers/bns_xcts/bns_xcts_regrid.hpp"
#include "Solvers/solver_startup.hpp"
#include <filesystem>
#include <type_traits>
namespace fs = std::filesystem;
using namespace Kadath;
using config_t = kadath_config_boost<BIN_INFO>;

// Foward declarations
int bns_xcts_driver (config_t& bconfig, std::string outputdir);
inline void bns_xcts_setup_space (config_t& bconfig);
inline void bns_xcts_setup_bin_config(config_t& bconfig);
inline void bns_xcts_superimposed_import(config_t& bconfig,std::array<std::string, 2> NSfilenames);
// end Forward declarations

int main(int argc, char** argv) {
  int rc = MPI_Init(&argc, &argv) ;
  if (rc!=MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl ;
    MPI_Abort(MPI_COMM_WORLD, rc) ;
  }
  int rank = 0 ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

  using InitSolver = Initialize_Solver<config_t>;
  
  // Initialize static member variables
  InitSolver::input_configname = "initial_bns.info";
  InitSolver::bconfig;
  InitSolver::rank = rank;

  // Run initialize routine based on CLI arguments
  InitSolver::init_solver(argc, argv);
  config_t bconfig = InitSolver::bconfig;
  
  if(InitSolver::example_setup) {
    if(rank == 0) {
      bconfig.initialize_binary({"ns","ns"});
      bconfig.set_defaults();
      bns_xcts_setup_bin_config(bconfig);
      bconfig.set_stage(TOTAL) = true;
      bconfig.write_config();
    }
  } else {
    bconfig.open_config();
    // solve at low res first    
    int const final_res = bconfig(BIN_RES);
    int const init_res = (std::isnan(bconfig.seq_setting(INIT_RES))) ? 
      9 : bconfig.seq_setting(INIT_RES);
    bool const res_inc = (init_res < final_res);
    auto const final_stages = bconfig.return_stages();
    auto & stages = bconfig.return_stages();
    
    if(res_inc) {
      // Solve only TOTAL before increasing resolution
      // All other stages only rescale the matter
      stages.fill(false);
      stages[TOTAL] = true;
    }
  
    if(InitSolver::setup_first) {
      bconfig.set_filename("initbin");
      bconfig.set(BIN_RES) = init_res;
      
      bns_xcts_setup_bin_config(bconfig);
      bns_xcts_setup_space(bconfig);
    }

    auto regrid = [&]() {
      std::string fname{"bns_regrid"};

      if(rank == 0)
        bns_xcts_regrid(bconfig, fname);
      bconfig.set_filename(fname);
      MPI_Barrier(MPI_COMM_WORLD);

      // Ensure input stages are activated
      stages = final_stages;
    };

    if(bconfig.control(REGRID) && !InitSolver::setup_first)
      regrid();
    
    int err = bns_xcts_driver(bconfig, InitSolver::outputdir);

    if(res_inc) {
      bconfig.set(BIN_RES) = final_res;
      regrid();
    
      // Rerun with new grid
      err = bns_xcts_driver(bconfig, InitSolver::outputdir);
    }
  }  
 
  MPI_Finalize();
  return EXIT_SUCCESS;
}
/**
 * bns_xcts_setup_bin
 *
 * Ensure proper binary settings
 *
 * @param[input] bconfig: BNS Configurator file
 */
inline void bns_xcts_setup_bin_config(config_t& bconfig) {
  check_dist(bconfig(DIST), bconfig(MADM, BCO1), bconfig(MADM, BCO2));
 
  // Binary Parameters
  bconfig.set(REXT) = 2 * bconfig(DIST);
  
  // equal mass system with "center of mass" at the origin
  bconfig.set(Q) = bconfig(MADM, BCO2) / bconfig(MADM, BCO1);
  
  // classical Newtonian estimate
  bconfig.set(COM) = 
    bco_utils::com_estimate(bconfig(DIST), bconfig(MADM, BCO1), bconfig(MADM, BCO2));
  
  // obtain 3PN estimate for the global, orbital omega
  bco_utils::KadathPNOrbitalParams(bconfig, \
        bconfig(MADM, BCO1), bconfig(MADM,BCO2));

  // delete ADOT, this can always be recalculated during
  // the eccentricity reduction stage
  bconfig.set(ADOT) = std::nan("1");
}

/**
 * bns_xcts_driver
 *
 * Control computation of BNS from setup config/dat file combination
 * filename is pulled from bconfig which should pair with <filename>.dat
 *
 * Also, KADATH at the time of writing, was strict on not allowing trivial
 * construction of base types (Base_tensor, Scalar, Tensor, Space, etc) which
 * means a driver is required to populate the related Solver class.
 *
 * @tparam[int] bconfig - Configurator file
 * @return Success/failure
 */
int bns_xcts_driver (config_t& bconfig, std::string outputdir){
  int exit_status = 0;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);

  std::string spacein = bconfig.space_filename();

  //just so you really know
  if(rank == 0) {
    std::cout << "Config File: " << bconfig.config_outputdir()+bconfig.config_filename() << std::endl
              << "Fields File: " << spacein << std::endl
              << bconfig  << std::endl;
  }
  FILE* ff1 = fopen (spacein.c_str(), "r") ;
  if(ff1 == NULL){
    std::cerr << spacein.c_str() << " failed to open for rank " << rank << "\n";
    std::_Exit(EXIT_FAILURE);
  }
  Space_bin_ns space (ff1) ;
	Scalar conf   (space, ff1) ;
	Scalar lapse  (space, ff1) ;
  Vector shift  (space, ff1) ;
  Scalar logh   (space, ff1) ;
  Scalar phi   (space, ff1) ;
	fclose(ff1) ;
  Base_tensor basis(space, CARTESIAN_BASIS);
  
  if(outputdir != "") {
    fs::create_directory(outputdir);
    bconfig.set_outputdir(outputdir) ;
  }
  if(bconfig.control(DELETE_SHIFT))
    shift.annule_hard();

  // load and setup the EOS
  const double h_cut = bconfig.template eos<double>(HCUT, BCO1);
  const std::string eos_file = bconfig.template eos<std::string>(EOSFILE, BCO1);
  const std::string eos_type = bconfig.template eos<std::string>(EOSTYPE, BCO1);

  if(eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
    bns_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
      bns_solver(bconfig, space, basis, conf, lapse, shift, logh, phi);
    bns_solver.solve();
  } else if(eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.template eos<int>(INTERP_PTS, BCO1) == 0) ? \
                            2000 : bconfig.template eos<int>(INTERP_PTS, BCO1);

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
    bns_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
      bns_solver(bconfig, space, basis, conf, lapse, shift, logh, phi);
    bns_solver.solve();
  } else { 
    std::cerr << "Unknown EOSTYPE." << endl;
    std::_Exit(EXIT_FAILURE);
  }

  return exit_status;
}

/**
 * bns_xcts_setup_space
 *
 * Create the numerical space and initial guess for the BNS by
 * obtaining the isolated TOV solutions, updating the BNS
 * config file, and importing the TOV solutions into the BNS space.
 *
 * @param[input] bconfig: BNS Configurator file
 */
inline void bns_xcts_setup_space (config_t& bconfig) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::array<int,2> bcos{BCO1, BCO2};
  std::array<std::string, 2> filenames;
  
  for(int i = 0; i < 2; ++i)
    filenames[i] = solve_NS_from_binary(bconfig, bcos[i]);
  
  // debugging only
  for(auto& f : filenames)
    if(rank == 0)
      std::cout << f << std::endl;

  if(rank == 0)
    bns_xcts_superimposed_import(bconfig, filenames);
  MPI_Barrier(MPI_COMM_WORLD);
  bconfig.control(SEQUENCES) = false;
  bconfig.open_config();
}

/**
 * bns_xcts_superimposed_import
 *
 * Flow control to initialize the binary space based on
 * superimposing two isolated compact object solutions.
 * These are simply read from file here before being passed
 * to bns_xcts_setup_boosted_3d for computing the initial
 * guess
 *
 * @param[input] bconfig: binary configurator
 * @param[input] NSfilenames: array of filenames for isolated solutions
 */
inline void bns_xcts_superimposed_import(config_t& bconfig,
  std::array<std::string, 2> NSfilenames) { 
  // load single NS configuration
  std::string ns1filename{NSfilenames[0]};
  kadath_config_boost<BCO_NS_INFO> NS1config(ns1filename);
  
  std::string ns2filename{NSfilenames[1]};
  kadath_config_boost<BCO_NS_INFO> NS2config(ns2filename);

  // setup eos and update central density
  const double h_cut = NS1config.eos<double>(HCUT);
  const std::string eos_file = NS1config.eos<std::string>(EOSFILE);
  const std::string eos_type = NS1config.eos<std::string>(EOSTYPE);

  if(eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
    bns_setup_boosted_3d<eos_t>(NS1config, NS2config, bconfig);
  } else if(eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (NS1config.eos<int>(INTERP_PTS) == 0) ? \
                            2000 : NS1config.eos<int>(INTERP_PTS);

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
    bns_setup_boosted_3d<eos_t>(NS1config, NS2config, bconfig);
  }
}
