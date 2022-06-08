/*
 * Copyright 2021
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
#include "kadath_bin_bh.hpp"
#include "mpi.h"
#include "Solvers/solver_startup.hpp"
#include "Solvers/bbh_xcts/bbh_xcts_solver.hpp"
#include "Solvers/bbh_xcts/bbh_xcts_regrid.hpp"
#include <filesystem>
#include <type_traits>
namespace fs = std::filesystem;
using namespace Kadath;
using config_t = kadath_config_boost<BIN_INFO>;

// Foward declarations
int bbh_xcts_driver (config_t& bconfig, std::string outputdir="");
inline void bbh_xcts_setup_space (config_t& bconfig);
inline void bbh_xcts_setup_bin_config(config_t& bconfig);
inline void bbh_xcts_superimposed_import(config_t& bconfig, std::array<std::string, 2> BHfilenames);
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
  InitSolver::input_configname = "initial_bbh.info";
  InitSolver::bconfig;
  InitSolver::rank = rank;

  // Run initialize routine based on CLI arguments
  InitSolver::init_solver(argc, argv);
  config_t bconfig = InitSolver::bconfig;

  if(InitSolver::example_setup) {
    if(rank == 0) {
      bconfig.initialize_binary({"bh","bh"});
      bconfig.set_defaults();
      bbh_xcts_setup_bin_config(bconfig);
      bconfig.write_config();
    }
  } else {
    bconfig.open_config();
    // solve at low res first    
    int const final_res = bconfig(BIN_RES);
    int const init_res = (std::isnan(bconfig.seq_setting(INIT_RES))) ? 
      9 : bconfig.seq_setting(INIT_RES);
    
    bool const res_inc = (final_res != init_res);
  
    if(InitSolver::setup_first) {
      bconfig.set_filename("initbin");
      bconfig.set(BIN_RES) = init_res;
      bbh_xcts_setup_bin_config(bconfig);
      bbh_xcts_setup_space(bconfig);
    }

    auto regrid = [&]() {
      std::string fname{"bbh_regrid"};

      if(rank == 0)
        bbh_xcts_regrid(bconfig, fname);
      bconfig.set_filename(fname);
      MPI_Barrier(MPI_COMM_WORLD);
    };
    
    if(bconfig.control(REGRID) && !InitSolver::setup_first)
      regrid();

    int err = bbh_xcts_driver(bconfig, InitSolver::outputdir);

    if(res_inc) {
      bconfig.set(BIN_RES) = final_res;
      regrid();
    
      // Rerun with new grid
      err = bbh_xcts_driver(bconfig, InitSolver::outputdir);
    }
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}

/**
 * bbh_xcts_setup_bin
 *
 * Ensure proper binary settings
 *
 * @param[input] bconfig: BBH Configurator file
 */
inline void bbh_xcts_setup_bin_config(config_t& bconfig){
  check_dist(bconfig(DIST), bconfig(MCH, BCO1), bconfig(MCH, BCO2));

  // Binary Parameters
  bconfig.set(REXT) = 2 * bconfig(DIST);
  bconfig.set(Q) = bconfig(MCH, BCO2) / bconfig(MCH, BCO1);
  
  // classical Newtonian estimate
  bconfig.set(COM) = 
    bco_utils::com_estimate(bconfig(DIST), bconfig(MCH, BCO1), bconfig(MCH, BCO2));
  
  // obtain 3PN estimate for the global, orbital omega
  bco_utils::KadathPNOrbitalParams(bconfig, \
        bconfig(MCH, BCO1), bconfig(MCH,BCO2));

  // delete ADOT, this can always be recalculated during
  // the eccentricity reduction stage
  bconfig.set(ADOT) = std::nan("1");
}

/**
 * bbh_xcts_driver
 *
 * Control computation of BBH from setup config/dat file combination
 * filename is pulled from bconfig which should pair with <filename>.dat
 *
 * Also, KADATH at the time of writing, was strict on not allowing trivial
 * construction of base types (Base_tensor, Scalar, Tensor, Space, etc) which
 * means a driver is required to populate the related Solver class.
 *
 * @return Success/failure
 */
int bbh_xcts_driver (config_t& bconfig, std::string outputdir) {
  int exit_status = 0;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // make sure outputdir directory exists for outputs
  if(rank == 0)
    std::cout << "Solutions will be stored in: " << outputdir << "\n" \
              << "Directory will be created if it doesn't exist.\n";

  std::string spacein = bconfig.space_filename();
  if(!fs::exists(spacein)) {
    // For debugging MPI bugs
    if(rank == 0) {
      std::cerr << "File: " << spacein << " not found.\n\n";
    } else {
      std::cerr << "File: " << spacein << " not found for another rank.\n\n";
    }
    std::_Exit(EXIT_FAILURE);
  }

  // just so you really know
  if(rank == 0) {
    std::cout << "Config File: " 
              << bconfig.config_outputdir()+bconfig.config_filename() << std::endl
              << "Fields File: " << spacein << std::endl
              << bconfig << std::endl;
  }
  FILE* ff1 = fopen (spacein.c_str(), "r") ;
  if(ff1 == NULL){
    // mainly for debugging MPI bugs
    std::cerr << spacein.c_str() << " failed to open for rank " << rank << "\n";
    std::_Exit(EXIT_FAILURE);
  }
  Space_bin_bh space (ff1) ;
	Scalar conf   (space, ff1) ;
	Scalar lapse  (space, ff1) ;
  Vector shift  (space, ff1) ;
	fclose(ff1) ;
  Base_tensor basis(space, CARTESIAN_BASIS);
  
  if(outputdir != "") {
    fs::create_directory(outputdir);
    bconfig.set_outputdir(outputdir) ;
  }
  if(bconfig.control(DELETE_SHIFT))
    shift.annule_hard();

  bbh_xcts_solver<decltype(bconfig), decltype(space)> 
      bbh_solver(bconfig, space, basis, conf, lapse, shift);
  bbh_solver.solve();
  
  return exit_status;
}

/**
 * bbh_xcts_setup_space
 *
 * Create the numerical space and initial guess for the BBH by
 * obtaining the isolated BH solutions, updating the BBH
 * config file, and importing the BH solutions into the BBH space.
 *
 * @param[input] bconfig: BBH Configurator file
 */
inline void bbh_xcts_setup_space (config_t& bconfig) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::array<int,2> bcos{BCO1, BCO2};
  std::array<std::string, 2> filenames;

  for(int i = 0; i < 2; ++i)
    filenames[i] = solve_BH_from_binary(bconfig, bcos[i]);  

  if(rank == 0)
    bbh_xcts_superimposed_import(bconfig, filenames);
  MPI_Barrier(MPI_COMM_WORLD);

  bconfig.open_config();
  bconfig.control(SEQUENCES) = false;
}

/**
 * bbh_xcts_superimposed_import
 *
 * Flow control to initialize the binary space based on
 * superimposing two isolated compact object solutions.
 * These are simply read from file here before being passed
 * to bbh_xcts_setup_boosted_3d for computing the initial
 * guess
 *
 * @param[input] bconfig: binary configurator
 * @param[input] BHfilenames: array of filenames for isolated solutions
 */
inline void bbh_xcts_superimposed_import(config_t& bconfig, std::array<std::string, 2> BHfilenames) {
  // load single BH configuration
  std::string bh1filename{BHfilenames[0]};
  kadath_config_boost<BCO_BH_INFO> BH1config(bh1filename);
  
  std::string bh2filename{BHfilenames[1]};
  kadath_config_boost<BCO_BH_INFO> BH2config(bh2filename);

  bbh_xcts_setup_boosted_3d(BH1config, BH2config, bconfig);
}