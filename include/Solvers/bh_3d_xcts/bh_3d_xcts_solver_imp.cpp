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
#include "bco_utilities.hpp"
#include "Solvers/co_solver_utils.hpp"
#include "bh_3d_xcts_solver.hpp"
#include "bh_3d_xcts_regrid.hpp"

using namespace Kadath;
using namespace Kadath::Margherita;

template<typename config_t, typename space_t>
bh_3d_xcts_solver<config_t, space_t>::bh_3d_xcts_solver(config_t& config_in, 
  space_t& space_in, Base_tensor& base_in,  
    Scalar& conf_in, Scalar& lapse_in, Vector& shift_in) :
      Solver<config_t, space_t>(config_in, space_in, base_in), 
        conf(conf_in), lapse(lapse_in), shift(shift_in), 
          fmet(Metric_flat(space_in, base_in))
{
  // initialize only the vector fields we need
  coord_vectors = default_co_vector_ary(space);

  update_fields_co(cfields, coord_vectors, {}, 0.);
}

// standardized filename for each converged dataset at the end of each stage.
template<typename config_t, typename space_t>
std::string bh_3d_xcts_solver<config_t, space_t>::converged_filename(
  const std::string stage) const {
  auto res = space.get_domain(0)->get_nbr_points()(0);
  std::stringstream ss;
  ss << "BH";
  if(stage != "") ss  << "_" << stage << ".";
  else ss << ".";
  ss << bconfig(MCH) << "." 
     << bconfig(CHI)<< "."
     << std::setfill('0')  << std::setw(1) << bconfig(NSHELLS) << "."
     << std::setfill('0')  << std::setw(2) << res;
  return ss.str();
}

template<typename config_t, typename space_t>
int bh_3d_xcts_solver<config_t, space_t>::solve() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int exit_status = EXIT_SUCCESS;
  
  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);

  //check if we have a non-boosted solution
  auto current_file = bconfig.config_filename_abs();
  auto const & final_chi = bconfig.seq_setting(FINAL_CHI);
  
  if(stage_enabled[TOTAL]) {
    exit_status = fixed_lapse_stage();
    if(exit_status == RELOAD_FILE)
      return exit_status;
  }

  if(stage_enabled[TOTAL_BC]) {    
    if(bconfig.control(ITERATIVE_CHI)) {      
      while(bconfig.control(ITERATIVE_CHI)) {        
        
        #ifdef DEBUG
        if(rank == 0) cout << bconfig << endl;
        #endif

        exit_status = von_Neumann_stage();
        if(exit_status == RELOAD_FILE)
          return exit_status;
        
        if((std::abs(final_chi) <= 0.8)
          || (std::abs(bconfig(CHI)) >= 0.8)) {
          bconfig.control(ITERATIVE_CHI) = false;
        } else if(std::abs(bconfig(CHI)) < 0.8)
          bconfig(CHI) = std::copysign(0.8, final_chi);
      }
      bconfig(CHI) = final_chi;
    }
    exit_status = von_Neumann_stage();
    if(exit_status == RELOAD_FILE)
      return exit_status;
  }
  
  // Barrier needed in case we need to read from the previous output
  MPI_Barrier(MPI_COMM_WORLD);
  if(last_stage_idx == BIN_BOOST)
    exit_status = RUN_BOOST;
  return exit_status;

}

template<typename config_t>
inline int bh_3d_xcts_driver (config_t& bconfig, std::string outputdir,
  kadath_config_boost<BIN_INFO> binconfig, const size_t bco){
  int exit_status = RELOAD_FILE;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(std::abs( bconfig(CHI)) > 0.84) {
    std::cerr << "Unable to handle chi > 0.84\n";
    return EXIT_FAILURE;
  }
  bconfig.seq_setting(INIT_RES) = (std::isnan(bconfig.seq_setting(INIT_RES))) ? 
    9 : bconfig.seq_setting(INIT_RES);

  // solve at low res first ? probably not needed
  const int final_res = bconfig(BCO_RES);
  bool res_inc = (bconfig.seq_setting(INIT_RES) < final_res);
  bconfig.set(BCO_RES) = bconfig.seq_setting(INIT_RES);

  // In the event we wish to solve for a highly spinning solution
  // we need to do an initial slow rotating solution before going to 
  // faster rotations otherwise the solution will diverge.
  bconfig.control(ITERATIVE_CHI) = std::abs(bconfig(CHI)) > 0.5;
  
  // We need to store the final desired CHI in case of iterative chi
  bconfig.seq_setting(FINAL_CHI) = bconfig(CHI);
  // Lower chi in case of iterative chi
  bconfig(CHI) = (bconfig.control(ITERATIVE_CHI)) ? 
    std::copysign(0.5, bconfig.seq_setting(FINAL_CHI)) : bconfig(CHI);
 
  // make sure BH directory exists for outputs
  if(rank == 0)
    std::cout << "Solutions will be stored in: " << outputdir << "\n" \
              << "Directory will be created if it doesn't exist.\n";
  fs::create_directory(outputdir);
  
  if(bconfig.control(SEQUENCES)) {
    bconfig.set_filename("initbh");
    if(rank == 0) {
      setup_co<BH>(bconfig);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // make sure all ranks have the same config
    bconfig.open_config();
    // FIXME - should this be turned off?
    bconfig.control(SEQUENCES) = false;
  }

  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);

  std::string spacein = bconfig.space_filename();
  if(!fs::exists(spacein)) {
    // mainly for debugging MPI bugs
    if(rank == 0) {
      std::cerr << "File: " << spacein << " not found.\n\n";
    } else {
      std::cerr << "File: " << spacein << " not found for another rank.\n\n";
    }
    std::_Exit(EXIT_FAILURE);
  }
  
  while(exit_status == RELOAD_FILE || exit_status == RUN_BOOST) { 
    spacein = bconfig.space_filename();
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
    Space_adapted_bh space (ff1) ;
    Scalar conf   (space, ff1) ;
    Scalar lapse  (space, ff1) ;
    Vector shift  (space, ff1) ;
    fclose(ff1) ;
    Base_tensor basis(space, CARTESIAN_BASIS);
    
    if(outputdir != "") bconfig.set_outputdir(outputdir) ;
    if(bconfig.control(DELETE_SHIFT))
      shift.annule_hard();

    bh_3d_xcts_solver<decltype(bconfig), decltype(space)> 
        bh_solver(bconfig, space, basis, conf, lapse, shift);
    exit_status = bh_solver.solve();
    if(exit_status == RUN_BOOST) {
      exit_status = bh_solver.binary_boost_stage(binconfig, bco);
    } else if(res_inc && exit_status != RELOAD_FILE) {
      res_inc = false;
      
      if(rank == 0)
        exit_status = bh_3d_xcts_regrid(bconfig, final_res, "initbh");
      bconfig.set_filename("initbh");
      bconfig.control(ITERATIVE_CHI) = false;
      MPI_Barrier(MPI_COMM_WORLD);
      bconfig.open_config();
      exit_status = RELOAD_FILE;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return exit_status;
}

template<typename config_t, typename space_t>
void bh_3d_xcts_solver<config_t, space_t>::syst_init(System_of_eqs& syst) {
  
  const int ndom = space.get_nbr_domains();
  // call the (flat) conformal metric "f"
  fmet.set_system(syst, "f");
 
  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_QPIG));
  syst.add_cst("PI"  , M_PI) ;
  syst.add_cst("M"   , bconfig(MIRR)) ;
  syst.add_cst("CM"  , bconfig(MCH));
  syst.add_cst("chi" , bconfig(CHI)) ;
  syst.add_var("ome" , bconfig(OMEGA));

  // include the coordinate fields
  syst.add_cst("mg"  , *coord_vectors[GLOBAL_ROT]);
  syst.add_cst("sm"  , *coord_vectors[S_BCO1]);
  
  syst.add_cst("ex"  , *coord_vectors[EX])  ;
  syst.add_cst("ey"  , *coord_vectors[EY])  ;
  syst.add_cst("einf", *coord_vectors[S_INF]);
  
  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("P"   , conf);
  syst.add_var("N"   , lapse);
  syst.add_var("bet" , shift) ;
  
  // define common combinations of conformal factor and lapse
  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");
  syst.add_def("A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");
  
  // definitions of integrals on the excision surface
  syst.add_def("intMsq= P^4 / 16. / PI") ;

  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
  syst.add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");
}

template<typename config_t, typename space_t>
void bh_3d_xcts_solver<config_t, space_t>::print_diagnostics_norot(const System_of_eqs & syst, 
    const int ite, const double conv) const {

	int ndom = space.get_nbr_domains() ;
  double r = bco_utils::get_radius(space.get_domain(1), OUTER_BC);
  
  Val_domain integMsq(syst.give_val_def("intMsq")()(2));
  double Mirrsq = space.get_domain(2)->integ(integMsq, INNER_BC);
  double Mirr = std::sqrt(Mirrsq);
  
  // compute the ADM mass as surface integral at infinity  
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  // compute the Komar mass as surface integral at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  // output to standard output  
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << "=======================================" << std::endl
            << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << std::endl
            << FORMAT << "Madm: " << Madm << std::endl
            << FORMAT << "Mk: " << Mk << " [" 
            << std::abs(Madm - Mk) / Madm << "]" << std::endl
            << FORMAT << "Mirr: " << Mirr << std::endl;
  std::cout << FORMAT << "R: " << r << std::endl;
  std::cout.flags(f);
} // end print diagnostics norot

// runtime diagnostics specific for rotating solutions
template<typename config_t, typename space_t>
void bh_3d_xcts_solver<config_t, space_t>::print_diagnostics(System_of_eqs const & syst, 
    const int ite, const double conv) const {

  // print all the diagnostics as in the non-rotating case first  
  print_diagnostics_norot(syst, ite, conv);

	int ndom = space.get_nbr_domains() ;
  
  Val_domain integS(syst.give_val_def("intS")()(2));
  double S = space.get_domain(2)->integ(integS, INNER_BC);

  Val_domain integMsq(syst.give_val_def("intMsq")()(2));
  double Mirrsq = space.get_domain(2)->integ(integMsq, INNER_BC);
  double Mch  = std::sqrt( Mirrsq + S * S / 4. / Mirrsq );

  // output the dimensionless spin and angular frequency parameter
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << FORMAT << "Omega: " << bconfig(OMEGA) << std::endl
            << FORMAT << "S: " << S << std::endl
            << FORMAT << "Mch: " << Mch << std::endl
            << FORMAT << "Chi: " << S / Mch / Mch << std::endl;
  std::cout.flags(f);
  std::cout << "=======================================" << "\n\n";
} // end print diagnostics rot
