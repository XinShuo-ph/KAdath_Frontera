#include "Solvers/solvers.hpp"
#include "mpi.h"
#include "bco_utilities.hpp"
#include "ns_3d_xcts_regrid.hpp"
#include <cmath>

using namespace Kadath;
using namespace Kadath::Margherita;

template<class eos_t, typename config_t, typename space_t>
ns_3d_xcts_solver<eos_t, config_t, space_t>::ns_3d_xcts_solver(config_t& config_in, 
  space_t& space_in, Base_tensor& base_in,  
    Scalar& conf_in, Scalar& lapse_in, Scalar& logh_in, Vector& shift_in) :
      Solver<config_t, space_t>(config_in, space_in, base_in), 
        conf(conf_in), lapse(lapse_in), logh(logh_in), shift(shift_in), 
          fmet(Metric_flat(space_in, base_in))
{
  // initialize only the vector fields we need
  coord_vectors = default_co_vector_ary(space);

  update_fields_co(cfields, coord_vectors,{}, 0.);
}

// standardized filename for each converged dataset at the end of each stage.
template<class eos_t, typename config_t, typename space_t>
std::string ns_3d_xcts_solver<eos_t, config_t, space_t>::converged_filename(
  const std::string& stage) const {
  auto res = space.get_domain(0)->get_nbr_points()(0);
  const std::string eos_file = bconfig.template eos<std::string>(EOSFILE);
  const std::string eosname = eos_file.substr(0, eos_file.find("."));
  std::stringstream ss;
  ss << "converged_NS";
  if(stage != "") ss  << "_" << stage << ".";
  else ss << ".";
  ss << eosname << "."
     << bconfig(MADM) << "."; 
  if(stage != "NOROT_BC") ss << bconfig(CHI)<< ".";
  else ss << "0.";
  ss << bconfig(NSHELLS) << "."
     <<std::setfill('0') << std::setw(2) << res;
  return ss.str();
}

template<class eos_t, typename config_t, typename space_t>
int ns_3d_xcts_solver<eos_t, config_t, space_t>::solve() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int exit_status = EXIT_SUCCESS;
  
  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);
  if(rank == 0) std::cout << "Last stage: " << last_stage << "\n";

  double const & final_chi = bconfig.seq_setting(FINAL_CHI);
  double const initial_chi = bconfig(CHI);
  
  if(stage_enabled[NOROT_BC]) {
    exit_status = norot_stage(false);
    bconfig(CHI) = initial_chi;
    return RELOAD_FILE; // FIXME: testing saving NOROT_BC ID for iterative boosts
  }

  if(stage_enabled[TOTAL_BC]) {
    bconfig.set_stage(NOROT_BC) = false;
    if(bconfig.control(ITERATIVE_CHI)) {
      exit_status = uniform_rot_stage();
      bconfig.control(ITERATIVE_CHI) = false;
      if(exit_status == RELOAD_FILE){        
        return exit_status;
      }
    }
    bconfig(CHI) = final_chi;
    exit_status = uniform_rot_stage();
    if(exit_status == RELOAD_FILE) {
      return exit_status;
    }
  }

  // Barrier needed in case we need to read from the previous output
  MPI_Barrier(MPI_COMM_WORLD);
  return exit_status;
}

template<typename config_t>
config_t ns_3d_xcts_isolated_driver (config_t& bconfig, std::string outputdir) {
  int exit_status = RELOAD_FILE;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  auto iterative_config = bconfig;

  // solve at low res first
  const int final_res = bconfig(BCO_RES);
  bool res_inc = (bconfig.seq_setting(INIT_RES) < final_res);
  bconfig.set(BCO_RES) = bconfig.seq_setting(INIT_RES);

  // in the event the user wants an MADM > MTOV
  const double final_MADM = bconfig(MADM);
  bconfig.control(ITERATIVE_M) = false;

  // in the event we need to iterate chi
  bconfig.seq_setting(FINAL_CHI) = bconfig(CHI);
  bconfig.control(ITERATIVE_CHI) = std::fabs(bconfig.seq_setting(FINAL_CHI)) > 0.2;

  // Lower chi in case of iterative chi
  bconfig(CHI) = (bconfig.control(ITERATIVE_CHI)) ? 
    std::copysign(0.1, bconfig.seq_setting(FINAL_CHI)) : bconfig(CHI);

  // make sure NS directory exists for outputs
  if(rank == 0)
    std::cout << "Solutions will be stored in: " << outputdir << "\n" \
              << "Directory will be created if it doesn't exist.\n";
  fs::create_directory(outputdir);

  if(std::isnan(bconfig.set(MADM)) && 
    (bconfig(MB) == 0. || std::isnan(bconfig.set(MB)))){
    if(rank == 0)
      std::cout << "Config error.  No madm nor mb found. \n\n";
    std::_Exit(EXIT_FAILURE);
  }
  
  if(bconfig.control(SEQUENCES)) {
    bconfig.set_filename("initns");
    if(rank == 0) {
      setup_co<NS>(bconfig);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // make sure all ranks have the same config
    bconfig.open_config();

    bconfig.control(ITERATIVE_M) = (bconfig(MADM) < final_MADM);
    if(bconfig.control(ITERATIVE_M) && std::fabs(bconfig(CHI)) < 1e-5) {
      if(rank == 0)
      std::cerr << "Cannot solve TOV for Madm = " << final_MADM
                << " without spin.\n";
      std::_Exit(EXIT_FAILURE);
    }
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

  // local version - need to save the control value in case of boost
  bool iter_m = bconfig.control(ITERATIVE_M);

  while(exit_status == RELOAD_FILE) {
    spacein = bconfig.space_filename();
    // just so you really know
    if(rank == 0) {
      std::cout << "Config File: " 
                << bconfig.config_filename_abs() << std::endl
                << "Fields File: " << spacein << std::endl
                << bconfig << std::endl;
    }
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    if(ff1 == NULL){
      // mainly for debugging MPI bugs
      std::cerr << spacein.c_str() << " failed to open for rank " << rank << "\n";
      std::_Exit(EXIT_FAILURE);
    }
    Space_spheric_adapted space (ff1) ;
    Scalar conf   (space, ff1) ;
    Scalar lapse  (space, ff1) ;
    Vector shift  (space, ff1) ;
    Scalar logh   (space, ff1) ;
    fclose(ff1) ;
    Base_tensor basis(space, CARTESIAN_BASIS);
    
    if(outputdir != "") bconfig.set_outputdir(outputdir) ;
    if(bconfig.control(DELETE_SHIFT))
      shift.annule_hard();

    // load and setup the EOS
    const double h_cut = bconfig.template eos<double>(HCUT);
    const std::string eos_file = bconfig.template eos<std::string>(EOSFILE);
    const std::string eos_type = bconfig.template eos<std::string>(EOSTYPE);

    if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
      ns_3d_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
        ns_solver(bconfig, space, basis, conf, lapse, logh, shift);
      exit_status = ns_solver.solve();

    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.template eos<int>(INTERP_PTS) == 0) ? \
                              2000 : bconfig.template eos<int>(INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      ns_3d_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
        ns_solver(bconfig, space, basis, conf, lapse, logh, shift);
      
      exit_status = ns_solver.solve();
    } else { 
      std::cerr << "Unknown EOSTYPE." << endl;
      std::_Exit(EXIT_FAILURE);
    }
    if(bconfig.set_stage(NOROT_BC) == true && exit_status == RELOAD_FILE){
      bconfig.set_stage(NOROT_BC) = false;
      iterative_config = bconfig; // save TOV solution before resolving
    } else if(bconfig.control(ITERATIVE_M) && exit_status != RELOAD_FILE) {
      
      bconfig.control(ITERATIVE_M) = false;
      bconfig.set(MADM) = final_MADM;
      bconfig.set_stage(NOROT_BC) = false;
      exit_status = RELOAD_FILE;
    } else if(res_inc && exit_status != RELOAD_FILE) {
      
      if(last_stage_idx != NOROT_BC)
        bconfig.set_stage(NOROT_BC) = false;
      
      bconfig(BCO_RES) += 2;
      res_inc = (bconfig(BCO_RES) < final_res);

      if(rank == 0)
        exit_status = ns_3d_xcts_interpolate_on_new_grid(bconfig, bconfig(BCO_RES), "initns");
      bconfig.set_filename("initns");
      MPI_Barrier(MPI_COMM_WORLD);
      bconfig.open_config();
      exit_status = RELOAD_FILE;
    } else if(bconfig.set_stage(NOROT_BC) == false && bconfig.set_stage(TOTAL_BC) == false){
      exit_status = EXIT_SUCCESS;
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
  // reset in case of boost stage
  bconfig.control(ITERATIVE_M) = iter_m;

  return iterative_config;
}

template<typename config_t>
inline int ns_3d_xcts_boosted_driver (config_t& nsconfig, config_t& iterative_config, 
  std::string outputdir, kadath_config_boost<BIN_INFO> binconfig, const size_t bco) {
  int exit_status = RUN_BOOST;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  auto bconfig = nsconfig;
  
  /* DNE
  auto bconfig = (nsconfig.control(ITERATIVE_M)) ? iterative_config : nsconfig;
  bool regridded = true;
  { // quick regrid and resolve
    if(rank == 0)
      cout << "Regridding before BOOST\n";
    std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
    auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);
    
    // ensure we don't redo norot for spinning data
    //if(last_stage_idx != NOROT_BC)
    bconfig.set_stage(NOROT_BC) = true;
    bconfig.set_stage(TOTAL_BC) = false;
    
    // regrid
    if(rank == 0)
      exit_status = ns_3d_xcts_interpolate_on_new_grid(bconfig, bconfig(BCO_RES), "initns");
    // update all ranks
    bconfig.set_filename("initns");
    MPI_Barrier(MPI_COMM_WORLD);
    bconfig.open_config();

    // controls to ensure the solution is solved based on new grid
    // without changing resolution
    bconfig.control(RESOLVE) = true;
    bconfig.seq_setting(INIT_RES) = bconfig(BCO_RES);
    bconfig.control(SEQUENCES) = false;
    ns_3d_xcts_isolated_driver (bconfig, outputdir);

    // reset controls
    bconfig.control(RESOLVE) = false;
    bconfig.seq_setting(INIT_RES) = 9;
    MPI_Barrier(MPI_COMM_WORLD);
    exit_status = RUN_BOOST;
    bconfig.set(CHI) = 0;
  }*/
  
  while(exit_status == RUN_BOOST) {
    auto spacein = bconfig.space_filename();

    // just so you really know
    if(rank == 0) {
      std::cout << "Config File: " 
                << bconfig.config_filename_abs() << std::endl
                << "Fields File: " << spacein << std::endl
                << bconfig << std::endl;
    }
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    if(ff1 == NULL){
      // mainly for debugging MPI bugs
      std::cerr << spacein.c_str() << " failed to open for rank " << rank << "\n";
      std::_Exit(EXIT_FAILURE);
    }
    Space_spheric_adapted space (ff1) ;
    Scalar conf   (space, ff1) ;
    Scalar lapse  (space, ff1) ;
    Vector shift  (space, ff1) ;
    Scalar logh   (space, ff1) ;
    fclose(ff1) ;
    Base_tensor basis(space, CARTESIAN_BASIS);
    
    if(outputdir != "") bconfig.set_outputdir(outputdir) ;
    if(bconfig.control(DELETE_SHIFT))
      shift.annule_hard();

    // load and setup the EOS
    const double h_cut = bconfig.template eos<double>(HCUT);
    const std::string eos_file = bconfig.template eos<std::string>(EOSFILE);
    const std::string eos_type = bconfig.template eos<std::string>(EOSTYPE);

    if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
      ns_3d_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
        ns_solver(bconfig, space, basis, conf, lapse, logh, shift);
        
      exit_status = ns_solver.binary_boost_stage(binconfig, bco);
    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.template eos<int>(INTERP_PTS) == 0) ? \
                              2000 : bconfig.template eos<int>(INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      ns_3d_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
        ns_solver(bconfig, space, basis, conf, lapse, logh, shift);
      
      exit_status = ns_solver.binary_boost_stage(binconfig, bco);
    } else { 
      std::cerr << "Unknown EOSTYPE." << endl;
      std::_Exit(EXIT_FAILURE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /* DNE
    if(nsconfig.control(ITERATIVE_M) && exit_status != RELOAD_FILE){
      nsconfig.control(ITERATIVE_M) = false;
      bconfig.set(MB) = nsconfig(MB);      
      bconfig.set(CHI) = nsconfig(CHI);
      bconfig.set(MADM) = nsconfig(MADM);
      exit_status = RELOAD_FILE;
    }*/
  }// end obtaining the rotating solution
  nsconfig = bconfig;
  return exit_status;
}
  

template<typename config_t>
inline int ns_3d_xcts_driver (config_t& bconfig, std::string outputdir, 
 kadath_config_boost<BIN_INFO> binconfig, const size_t bco) {
  if(std::isnan(bconfig.seq_setting(INIT_RES))) 
    bconfig.seq_setting(INIT_RES) = 9;
  
  int exit_status = EXIT_SUCCESS;
  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);
  
  auto iterative_config = ns_3d_xcts_isolated_driver(bconfig, outputdir);

  if(last_stage_idx == BIN_BOOST)
    exit_status = ns_3d_xcts_boosted_driver(bconfig, iterative_config, outputdir, binconfig, bco);  

  return exit_status;
}

template<class eos_t, typename config_t, typename space_t>
void ns_3d_xcts_solver<eos_t, config_t, space_t>::syst_init(System_of_eqs& syst) {
  using namespace Kadath::Margherita;
  
  auto& space = syst.get_space();
  const int ndom = space.get_nbr_domains();
  // call the (flat) conformal metric "f"
  fmet.set_system(syst, "f");
 
  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_QPIG));
  
  // include the coordinate fields
  syst.add_cst("mg"  , *coord_vectors[GLOBAL_ROT]);
  syst.add_cst("sm"  , *coord_vectors[S_BCO1]);
  syst.add_cst("einf", *coord_vectors[S_INF]);
  
  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("P"   , conf);
  syst.add_var("N"   , lapse);
  
  // define common combinations of conformal factor and lapse
  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");
 
  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
  syst.add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");
  syst.add_def(ndom - 1, "intMadmalt = -dr(P) * 2 / 4piG");
  
  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");
 
  // define the EOS operators
  Param p;
  syst.add_ope("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope("rho", &EOS<eos_t,DENSITY>::action, &p);
  syst.add_ope("dHdlnrho", &EOS<eos_t,DHDRHO>::action, &p);
 
  // define rest-mass density, internal energy and pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");
  syst.add_def("dHdlnrho = dHdlnrho(h)");

  // definition to rescale the equations
  // delta = p / rho
  syst.add_def("delta = h - eps - 1.");

}

template<class eos_t, typename config_t, typename space_t>
void ns_3d_xcts_solver<eos_t, config_t, space_t>::print_diagnostics_norot(const System_of_eqs & syst, 
    const int ite, const double conv) const {

  // total number of domains	
  int ndom = space.get_nbr_domains() ;

  // compute the baryonic mass at volume integral from the given integrant
  double baryonic_mass =
      syst.give_val_def("intMb")()(0).integ_volume() +
      syst.give_val_def("intMb")()(1).integ_volume();

  // compute the ADM mass as surface integral at infinity  
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  // compute the Komar mass as surface integral at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  // get the maximum and minimum coordinate radius along the surface,
  // i.e. the adapted domain boundary
  auto rs = bco_utils::get_rmin_rmax(space, 1);

  // alternative, equivalent ADM mass integral
  Val_domain integMadmalt(syst.give_val_def("intMadmalt")()(ndom - 1));
  double Madmalt = space.get_domain(ndom - 1)->integ(integMadmalt, OUTER_BC);

  // output to standard output  
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << "=======================================" << std::endl
            << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << std::endl
            << FORMAT << "Mb: " << baryonic_mass << std::endl
            << FORMAT << "Madm: " << Madm << std::endl
            << FORMAT << "Madm_ql: " << Madmalt 
            << " [" << std::abs(Madm - Madmalt) / Madm << "]" << std::endl
            << FORMAT << "Mk: " << Mk << " [" 
            << std::abs(Madm - Mk) / Madm << "]" << std::endl;
  std::cout << FORMAT << "R: " << rs[0] << " " << rs[1] << "\n";
  std::cout.flags(f);
} // end print diagnostics norot

// runtime diagnostics specific for rotating solutions
template<class eos_t, typename config_t, typename space_t>
void ns_3d_xcts_solver<eos_t, config_t, space_t>::print_diagnostics(System_of_eqs const & syst, 
    const int ite, const double conv) const {

  // print all the diagnostics as in the non-rotating case first  
  print_diagnostics_norot(syst, ite, conv);

  // total number of domains
	int ndom = space.get_nbr_domains() ;

  // compute the ADM angular momentum as surface integral at infinity  
  Val_domain integJ(syst.give_val_def("intJ")()(ndom - 1));
  double J = space.get_domain(ndom - 1)->integ(integJ, OUTER_BC);

  // compute the ADM mass as surface integral at infinity    
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  // output the dimensionless spin and angular frequency parameter
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << FORMAT << "Jadm: " << J << std::endl
            << FORMAT << "Chi: " << J / Madm / Madm << " [" << bconfig(CHI) << "]\n"
            << FORMAT << "Omega: " << bconfig(OMEGA) << std::endl;
  std::cout.flags(f);
  std::cout << "=======================================" << "\n\n";
} // end print diagnostics rot

template<class eos_t, typename config_t, typename space_t>
void ns_3d_xcts_solver<eos_t, config_t, space_t>::update_config_quantities(const double& loghc) {
  bconfig.set(HC) = std::exp(loghc);
  bconfig.set(NC) = EOS<eos_t,DENSITY>::get(bconfig(HC));
}
