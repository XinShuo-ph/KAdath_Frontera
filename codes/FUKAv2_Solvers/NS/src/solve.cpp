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
#include "mpi.h"
#include "Solvers/ns_3d_xcts/ns_3d_xcts_solver.hpp"
#include "Solvers/solver_startup.hpp"
using namespace Kadath;

int main(int argc, char** argv) {
  int rc = MPI_Init(&argc, &argv) ;
  if (rc!=MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl ;
    MPI_Abort(MPI_COMM_WORLD, rc) ;
  }
  int rank = 0 ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
  using config_t = kadath_config_boost<BCO_NS_INFO>;
  using InitSolver = Initialize_Solver<config_t>;
  
  // Initialize static member variables
  InitSolver::input_configname = "initial_ns.info";
  InitSolver::bconfig;
  InitSolver::rank = rank;

  // Run initialize routine based on CLI arguments
  InitSolver::init_solver(argc, argv);
  config_t bconfig = InitSolver::bconfig;

  if(InitSolver::example_setup) {
    if(rank == 0) {
      // Generate <example name>.info and terminate      
      bconfig.set_defaults();
      bconfig.control(SEQUENCES) = InitSolver::setup_first;
      bconfig.write_config();
    }
  } else {
    bconfig.open_config();
    bconfig.control(SEQUENCES) = InitSolver::setup_first;
    int err = ns_3d_xcts_driver(bconfig, InitSolver::outputdir);
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}
