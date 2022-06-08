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
#include "Solvers/bco_solver_utils.hpp"
#include "Solvers/bh_3d_xcts/bh_3d_xcts_solver.hpp"
#include "Solvers/ns_3d_xcts/ns_3d_xcts_solver.hpp"
#include "kadath.hpp"

using namespace Kadath;
using namespace Kadath::Margherita;

template<class eos_t, typename config_t, typename space_t>
bhns_xcts_solver<eos_t, config_t, space_t>::bhns_xcts_solver(config_t& config_in, 
  space_t& space_in, Base_tensor& base_in,  
    Scalar& conf_in, Scalar& lapse_in, Vector& shift_in, Scalar& logh_in, Scalar& phi_in) :
      Solver<config_t, space_t>(config_in, space_in, base_in), 
        conf(conf_in), lapse(lapse_in), shift(shift_in), logh(logh_in), phi(phi_in),
          fmet(Metric_flat(space_in, base_in)),
            xc1(bco_utils::get_center(space_in,space.NS)),
              xc2(bco_utils::get_center(space_in,space.BH)),
                xo(bco_utils::get_center(space,ndom-1)),
                  excluded_doms({space.BH,space.BH+1})
{
  // initialize coordinate vector fields
  coord_vectors = default_binary_vector_ary(space);

  update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);
}

// standardized filename for each converged dataset at the end of each stage.
template<class eos_t, typename config_t, typename space_t>
std::string bhns_xcts_solver<eos_t, config_t, space_t>::converged_filename(
  const std::string& stage) const {
  const std::string eos_file = bconfig.template eos<std::string>(EOSFILE, BCO1);
  const std::string eosname = eos_file.substr(0, eos_file.find("."));
  std::stringstream ss;
  auto res = space.get_domain(0)->get_nbr_points()(0);
  auto M1 = bconfig(MADM, BCO1);
  auto M2 = bconfig(MCH, BCO2);
  if(M2 > M1) std::swap(M1, M2);
  bconfig.set(Q) = M2 / M1;
  auto Mtot = M1 + M2;
  ss << "converged_BHNS";
  if(stage != "") ss  << "_" << stage << ".";
  else ss << ".";
  ss << eosname << "."
     << bconfig(DIST)      << "." 
     << bconfig(CHI, BCO1) << "." 
     << bconfig(CHI, BCO2) << "."
     << Mtot << ".q" 
     << bconfig(Q)         << "."
     << std::setfill('0')  << std::setw(1) << bconfig(NSHELLS, BCO1) << "."
     << std::setfill('0')  << std::setw(1) << bconfig(NSHELLS, BCO2) << "."
     << std::setfill('0')  << std::setw(2) << res;
  return ss.str();
}

template<class eos_t, typename config_t, typename space_t>
int bhns_xcts_solver<eos_t, config_t, space_t>::solve() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int exit_status = EXIT_SUCCESS;

  auto stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);
  
  //check if we have solved this before
  if(solution_exists(last_stage)) {
    if(rank == 0)
      std::cout << "Solved previously: " << bconfig.config_filename_abs() << std::endl;
    return EXIT_SUCCESS;
  }
  
  // overview output
	if (rank==0) {
		std::cout << "=================================" << endl;
		std::cout << "BHNS grav input" << endl;
		std::cout << "Distance: " << bconfig(DIST) << endl;
		std::cout << "Omega guess: " << bconfig(GOMEGA) << endl;
		std::cout << "Units: " << bconfig(QPIG) << endl;
		std::cout << "=================================" << endl;
	}
  
  if(stage_enabled[TOTAL]) {
    exit_status = hydrostatic_equilibrium_stage(TOTAL, "TOTAL");
    if(exit_status != EXIT_SUCCESS) std::_Exit(EXIT_FAILURE);
  }
  
  if(stage_enabled[TOTAL_BC]) {
    if(bconfig.control(FIXED_GOMEGA)) {
      exit_status = hydro_rescaling_stages(TOTAL_BC, "TOTAL_BC_FIXED_OMEGA");
      
      // In the event we are creating a new BHNS sequence, FIXED_OMEGA
      // is used only for the initial solution before resolving the hydro
      // consistently.  The constrols are then disabled since we only need to
      // do this once
      if(exit_status == EXIT_SUCCESS && bconfig.control(SEQUENCES)) {
        bconfig.control(SEQUENCES) = false;
        bconfig.control(FIXED_GOMEGA) = false;

        exit_status = hydrostatic_equilibrium_stage(TOTAL_BC, "TOTAL_BC");
      }
    } else
      exit_status = hydrostatic_equilibrium_stage(TOTAL_BC, "TOTAL_BC");
    
    if(exit_status != EXIT_SUCCESS) std::_Exit(EXIT_FAILURE);
  }
  if(stage_enabled[ECC_RED]) {
    exit_status = hydro_rescaling_stages(ECC_RED, "ECC_RED");
    if(exit_status != EXIT_SUCCESS) std::_Exit(EXIT_FAILURE);
  }  
 
  // Barrier needed in case we need to read from the previous output
  MPI_Barrier(MPI_COMM_WORLD);
 
  return exit_status;

}

template<class eos_t, typename config_t, typename space_t>
void bhns_xcts_solver<eos_t, config_t, space_t>::syst_init(System_of_eqs& syst) {
  using namespace Kadath::Margherita;
  
  const int ndom = space.get_nbr_domains();
  // call the (flat) conformal metric "f"
  fmet.set_system(syst, "f");
  
  // define the EOS operators
  Param p;
  syst.add_ope("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope("rho", &EOS<eos_t,DENSITY>::action, &p);
  syst.add_ope("dHdlnrho", &EOS<eos_t,DHDRHO>::action, &p);
 
  // define numerical constants
  syst.add_cst("4piG", bconfig(QPIG));
  
  // baryonic mass and dimensionless spin are fixed input parameters
  // along with the ADM of each NS at infinite separation
  syst.add_cst("Madm1" , bconfig(MADM  , BCO1));
  syst.add_cst("Mb1" , bconfig(MB , BCO1));
  syst.add_cst("chi1", bconfig(CHI, BCO1));

  syst.add_cst("Mch" , bconfig(MCH , BCO2));
  syst.add_cst("Mirr", bconfig(MIRR, BCO2));
  syst.add_cst("chi2", bconfig(CHI, BCO2));

  // coordinate dependent (rotational) vector fields
  syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
  syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
  syst.add_cst("mp", *coord_vectors[BCO2_ROT]);

  // cartesianbasis vector fields
  syst.add_cst("ex", *coord_vectors[EX]);
  syst.add_cst("ey", *coord_vectors[EY]);
  syst.add_cst("ez", *coord_vectors[EZ]);

  // flat spac surface normals
  syst.add_cst("sm"  , *coord_vectors[S_BCO1]);
  syst.add_cst("sp"  , *coord_vectors[S_BCO2]);
  syst.add_cst("einf", *coord_vectors[S_INF]);
  
  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("P", conf);
  syst.add_var("N", lapse);
  syst.add_var("bet", shift);

  // quasi-local measurement of the component ADM masses
  syst.add_var("qlMadm1" , bconfig(QLMADM, BCO1));
  
  // check if corotation is considered and adjust
  // rotational velocity components accordingly
  if(bconfig.control(COROT_BIN)) {
    // for a corotating binary, there is no angular frequency parameter
    // since stellar rotation is fixed by the orbital frequency
    bconfig.set(OMEGA, BCO1) = 0.;
    syst.add_cst("omes1", bconfig(OMEGA,BCO1));

    bconfig.set(OMEGA, BCO2) = 0.;
    syst.add_cst("omes2", bconfig(OMEGA,BCO2));
  } else {
    // in the arbitrary spinning / irrotational case
    // the irrotational part is defined by the velocity potential
    syst.add_var("phi" , phi);

    // check if stellar omega is fixed or has to be solved
    // for a specific dimensionless spin
    if(!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
      syst.add_cst("omes1", bconfig(FIXED_BCOMEGA,BCO1));
      bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA,BCO1);
    } else {
      syst.add_var("omes1", bconfig(OMEGA,BCO1));
    }

    if(!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
      syst.add_cst("omes2", bconfig(FIXED_BCOMEGA,BCO2));
      bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA,BCO2);
    } else {
      syst.add_var("omes2", bconfig(OMEGA,BCO2));
    }

    // for mixed spins, we have additional definitions
    // that are required for both objects that describe
    // the local rotation field
    for(int d = space.NS; d <= space.ADAPTEDNS; ++d){
      syst.add_def(d, "s^i  = omes1 * mm^i");
      syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
    }
    syst.add_def(space.ADAPTEDBH+1, "s^i = omes2 * mp^i"); 
  }
  
  // define common combinations of conformal factor and lapse
  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");
 
  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
  syst.add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");
  syst.add_def(ndom - 1, "intMadmalt = -dr(P) * 2 / 4piG");
  syst.add_def("intMsq  = P^4 / 4. / 4piG") ;

  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");
  
  // define rest-mass density, internal energy and pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");
  syst.add_def("dHdlnrho = dHdlnrho(h)");
  
  // definition to rescale the equations
  // delta = p / rho
  syst.add_def("delta = h - eps - 1.");
 
  // the conformal extrinsic curvature
  syst.add_def("A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");
  
  // ADM linear momentum, surface integrant at infinity
  syst.add_def(ndom-1, "intPx = A_i^j * ex_j * einf^i");
  syst.add_def(ndom-1, "intPy = A_i^j * ey_j * einf^i");
  syst.add_def(ndom-1, "intPz = A_i^j * ez_j * einf^i");

  // quasi-local spin, surface integral outside the stellar matter distribution
  syst.add_def(space.ADAPTEDNS+1, "intS1 = A_ij * mm^i * sm^j / 2. / 4piG");
  syst.add_def(space.ADAPTEDBH+1, "intS2 = A_ij * mp^i * sp^j / 2. / 4piG");
}

// runtime diagnostics specific for rotating solutions
template<class eos_t, typename config_t, typename space_t>
void bhns_xcts_solver<eos_t, config_t, space_t>::print_diagnostics(System_of_eqs const & syst, 
    const int ite, const double conv) const {
  std::ios_base::fmtflags f( std::cout.flags() );
  
  double baryonic_massNS = 0.;
  double ql_massNS = 0.;
  for(int d = space.NS; d <= space.ADAPTEDNS; ++d){
    baryonic_massNS += syst.give_val_def("intMb")()(d).integ_volume();
    ql_massNS += syst.give_val_def("intM")()(d).integ_volume();
  }
  
  auto [ rmin, rmax ] = bco_utils::get_rmin_rmax(space, space.ADAPTEDNS);
  double r_bh = bco_utils::get_radius(space.get_domain(space.ADAPTEDBH), EQUI) ;
  
  double mirrsq = space.get_domain(space.BH+2)->integ(syst.give_val_def("intMsq")()(space.BH+2), INNER_BC);
  double Mirr = std::sqrt(mirrsq);

  Val_domain integS1(syst.give_val_def("intS1")()(space.ADAPTEDNS+1));
  double S1 = space.get_domain(space.ADAPTEDNS+1)->integ(integS1, OUTER_BC);
  double chi1 = S1 / bconfig(MADM, BCO1) / bconfig(MADM, BCO1);
  double chiql1 = S1 / ql_massNS / ql_massNS;

  Val_domain integS2(syst.give_val_def("intS2")()(space.ADAPTEDBH+1));
  double S2 = space.get_domain(space.ADAPTEDBH+1)->integ(integS2, OUTER_BC);
  double Mch  = std::sqrt( mirrsq + S2 * S2 / 4. / mirrsq );
  double chi2 = S2 / Mch / Mch;

  Val_domain integPx(syst.give_val_def("intPx")()(ndom - 1));
  double Px = space.get_domain(ndom - 1)->integ(integPx, OUTER_BC);

  Val_domain integPy(syst.give_val_def("intPy")()(ndom - 1));
  double Py = space.get_domain(ndom - 1)->integ(integPy, OUTER_BC);

  Val_domain integPz(syst.give_val_def("intPz")()(ndom - 1));
  double Pz = space.get_domain(ndom - 1)->integ(integPz, OUTER_BC);

  std::cout << "=======================================" << endl;
  std::cout << FORMAT << "Iter: "  << ite << std::endl
            << FORMAT << "Error: " << conv << "\n" ;
  std::cout << FORMAT << "Omega: " << bconfig(GOMEGA) << endl ;
  std::cout << FORMAT << "Axis: "  << bconfig(COM) << endl ;
  std::cout << FORMAT << "Padm: "  << "[" <<Px << ", " << Py << ", " << Pz << "]\n\n";
  
  std::cout << FORMAT << "NS-Mb: "    << baryonic_massNS << "[" << bconfig(MB,BCO1) << "]\n";
  std::cout << FORMAT << "NS-Madm_ql: " << ql_massNS << "[" << bconfig(MADM,BCO1) << "]\n";
  std::cout << FORMAT << "NS-R: "     << rmin   << " " << rmax << std::endl
            << FORMAT << "NS-S: "     << S1     << std::endl;
  std::cout << FORMAT << "NS-Chi: "   << chi1   << "[" << bconfig(CHI,BCO1) << "]\n";
  std::cout << FORMAT << "NS-qlChi: " << chiql1 << std::endl;
  std::cout << FORMAT << "NS-Omega: " << bconfig(OMEGA, BCO1) << "\n\n";

  std::cout << FORMAT << "BH-Mirr: "  << Mirr << "[" << bconfig(MIRR,BCO2) << "]\n";
  std::cout << FORMAT << "BH-Mch: "   << Mch  << "[" << bconfig(MCH,BCO2) << "]\n";
  std::cout << FORMAT << "BH-R: "     << r_bh << std::endl;
  std::cout << FORMAT << "BH-S: "     << S2   << std::endl;
  std::cout << FORMAT << "BH-Chi: "   << chi2 << "[" << bconfig(CHI,BCO2) << "]\n";
  std::cout << FORMAT << "BH-Omega: " << bconfig(OMEGA, BCO2) << std::endl;
  std::cout.flags(f);
  std::cout << "=======================================" << endl;
} // end print_diagnostics

/*template<class eos_t, typename config_t, typename space_t>
void ns_3d_xcts_solver<eos_t, config_t, space_t>::update_config_quantities(const double& loghc) {
  bconfig.set(HC) = std::exp(loghc);
  bconfig.set(NC) = EOS<eos_t,DENSITY>::get(bconfig(HC));
}*/

template<typename eos_t>
inline void bhns_setup_boosted_3d(kadath_config_boost<BCO_NS_INFO>& NSconfig, 
  kadath_config_boost<BCO_BH_INFO>& BHconfig,
  kadath_config_boost<BIN_INFO>& bconfig) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // open previous ns solution
  std::string nsspaceinf = NSconfig.space_filename();
  FILE* ff1 = fopen(nsspaceinf.c_str(), "r") ;
	Space_spheric_adapted nsspacein(ff1) ;
	Scalar nsconf   (nsspacein, ff1) ;
	Scalar nslapse  (nsspacein, ff1) ;
	Vector nsshift  (nsspacein, ff1) ;
  Scalar nslogh   (nsspacein, ff1) ;
  Scalar nsphi    (nsspacein, ff1) ;
	fclose(ff1) ;
  // end opening bns solution

  // update NSconfig quantities before updating binary configuration file
  NSconfig.set(HC) = std::exp(bco_utils::get_boundary_val(0, nslogh));
  NSconfig.set(NC) = EOS<eos_t,DENSITY>::get(NSconfig(HC)) ;
  bco_utils::update_config_NS_radii(nsspacein, NSconfig, 1);
  
  // obtain adapted NS shells for radius information and copying adapted mapping later
  const Domain_shell_outer_adapted* old_outer_adapted1 =
        dynamic_cast<const Domain_shell_outer_adapted*>(nsspacein.get_domain(1));
  const Domain_shell_inner_adapted* old_inner_adapted1 =
        dynamic_cast<const Domain_shell_inner_adapted*>(nsspacein.get_domain(2));
  
  // setup radius field - needed for copying the adapted domain mappings.
  Scalar old_space_radius(nsspacein);
  old_space_radius.annule_hard();

	int ndominns = nsspacein.get_nbr_domains() ;

  for(int d = 0; d < ndominns; ++d) 
    old_space_radius.set_domain(d) = nsspacein.get_domain(d)->get_radius();
  
  old_space_radius.set_domain(1) = old_outer_adapted1->get_outer_radius();
  old_space_radius.std_base();
  // end setup radius field

  // open old BH solution
  std::string bhspaceinf = BHconfig.space_filename(); 
  FILE* ff2 = fopen(bhspaceinf.c_str(), "r") ;
	Space_adapted_bh bhspacein(ff2) ;
	Scalar bhconf  (bhspacein, ff2) ;
	Scalar bhlapse (bhspacein, ff2) ;
	Vector bhshift (bhspacein, ff2) ;
	fclose(ff2) ;
  // end open BH solution
  bco_utils::update_config_BH_radii(bhspacein, BHconfig, 1, bhconf);
  
  auto interp_field = [&](auto& space, int outer_dom, auto& old_phi) {
    const int d = outer_dom;
    const Domain_shell_inner_adapted* old_inner = 
      dynamic_cast<const Domain_shell_inner_adapted*>(space.get_domain(d+1));
    //interpolate old_phi field outside of the star for import
    bco_utils::update_adapted_field(old_phi, d, d+1, old_inner, INNER_BC);
  };

  // in case we used boosted TOVs, we need to import PHI
  if(NSconfig.set_field(PHI) == true)
    interp_field(nsspacein, 1, nsphi);

  //start Update config vars
  double r_max_tot = std::max(bconfig(RMID, BCO1), bconfig(RMID,BCO2));
  const double rout_sep_est = (bconfig(DIST) / 2. - r_max_tot) / 3. + r_max_tot;
  const double rout_max_est = 1.5 * r_max_tot;
  bconfig.set(ROUT, BCO1) = (rout_sep_est > rout_max_est) ? rout_max_est : rout_sep_est;
  bconfig.set(ROUT, BCO2) = bconfig(ROUT,BCO1);
  //end updating config vars

  // setup domain boundaries
  std::vector<double> out_bounds(1+bconfig(OUTER_SHELLS));
  std::vector<double> NS_bounds(3+bconfig(NSHELLS,BCO1)+bconfig(NINSHELLS,BCO1));
  std::vector<double> BH_bounds(3+bconfig(NSHELLS,BCO2));
  
  bco_utils::set_BH_bounds(BH_bounds, bconfig, BCO2, true);
  bco_utils::set_NS_bounds(NS_bounds, bconfig, BCO1);

  //for out_bounds.size > 1 - add equi-distance shells
  for(int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(REXT) + e * 0.25 * bconfig(REXT);

  bco_utils::print_bounds("NS-bounds", NS_bounds);
  bco_utils::print_bounds("BH-bounds", BH_bounds);
  bco_utils::print_bounds("outer-bounds", out_bounds);

  //Setup actual space
  int typer = CHEB_TYPE ;
  Space_bhns space (typer, bconfig(DIST), NS_bounds, BH_bounds, out_bounds, bconfig(BIN_RES), bconfig(NINSHELLS, BCO1));
  Base_tensor basis (space, CARTESIAN_BASIS) ;

  const Domain_shell_inner_adapted*    new_ns_inner = dynamic_cast<const Domain_shell_inner_adapted*>(space.get_domain(space.ADAPTEDNS+1));
  const Domain_shell_outer_adapted*    new_ns_outer = dynamic_cast<const Domain_shell_outer_adapted*>(space.get_domain(space.ADAPTEDNS));
  
  const Domain_shell_outer_homothetic* old_bh_outer = dynamic_cast<const Domain_shell_outer_homothetic*>(bhspacein.get_domain(1));
  
  //Update BH fields based to help with interpolation later
  bco_utils::update_adapted_field(bhconf, 2, 1, old_bh_outer, OUTER_BC);
  bco_utils::update_adapted_field(bhlapse, 2, 1, old_bh_outer, OUTER_BC);
  
  //Updated mapping for NS
  bco_utils::interp_adapted_mapping(new_ns_inner, 1, old_space_radius);
  bco_utils::interp_adapted_mapping(new_ns_outer, 1, old_space_radius);
  
  double xc1 = bco_utils::get_center(space, space.NS);
  double xc2 = bco_utils::get_center(space, space.BH);

  std::cout << "xc1: " << xc1 << std::endl;
  std::cout << "xc2: " << xc2 << std::endl;
  
  if(NSconfig.set_field(PHI) == true)
    bco_utils::update_adapted_field(nsphi, 1, 2, old_inner_adapted1, INNER_BC);
  
  Scalar logh(space);
  logh.annule_hard();
	logh.std_base();
  
  Scalar conf(space);
  conf.annule_hard();

  Scalar lapse(space);
  lapse.annule_hard();

  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();

  Scalar phi(space);
  phi.annule_hard();

  const double ns_invw4 = bco_utils::set_decay(bconfig, BCO1);
  const double bh_invw4 = bco_utils::set_decay(bconfig, BCO2);
  
  if(rank == 0)
    std::cout << "WeightNS: " << bconfig(DECAY, BCO1) << ", "
              << "WeightBH: " << bconfig(DECAY, BCO2) << std::endl;
  
  int ndom = space.get_nbr_domains();
  for(int dom = 0; dom < ndom; dom++)
  {
    // get an index in each domain to iterate over all colocation points
		Index new_pos(space.get_domain(dom)->get_nbr_points());

		do {
      // get cartesian coordinates of the current colocation point
			double x = space.get_domain(dom)->get_cart(1)(new_pos);
			double y = space.get_domain(dom)->get_cart(2)(new_pos);
			double z = space.get_domain(dom)->get_cart(3)(new_pos);

      // define a point shifted suitably to the stellar centers in the binary
			Point absol1(3);
			absol1.set(1) = (x - xc1);
			absol1.set(2) = y;
			absol1.set(3) = z;
      double r2 = y * y + z * z; 
      double r2_1 = (x - xc1) * (x - xc1) + r2;
      double r4_1  = r2_1 * r2_1;
      double r4_invw4_1 = r4_1 * ns_invw4;
      double decay_1 = std::exp(-r4_invw4_1);

      Point absol2(3);
			absol2.set(1) = (x - xc2);
			absol2.set(2) = y;
			absol2.set(3) = z;
      double r2_2 = (x - xc2) * (x - xc2) + r2;
      double r4_2  = r2_2 * r2_2;
      double r4_invw4_2 = r4_2 * bh_invw4;
      double decay_2 = std::exp(-r4_invw4_2);
      
      if (dom < ndom - 1) {
        conf .set_domain(dom).set(new_pos) =              \
          1. + decay_1 * (nsconf.val_point(absol1) - 1.) \
             + decay_2 * (bhconf.val_point(absol2) - 1.);
        
        lapse.set_domain(dom).set(new_pos) =               \
          1. + decay_1 * (nslapse.val_point(absol1) - 1.) \
             + decay_2 * (bhlapse.val_point(absol2) - 1.);
        
        logh .set_domain(dom).set(new_pos) = 0. + \
          decay_1 * nslogh.val_point(absol1);
        
        phi.set_domain(dom).set(new_pos) = 0;
        if(NSconfig.set_field(PHI) == true)
          phi.set_domain(dom).set(new_pos) += \
            decay_1 * nsphi.val_point(absol1);

        for (int i = 1; i <= 3; i++)
          shift.set(i).set_domain(dom).set(new_pos) =   \
            0. + decay_1 * nsshift(i).val_point(absol1) \
               + decay_2 * bhshift(i).val_point(absol2);
   
      } else {
        // We have to set the compactified domain manually 
        // since the outer collocation point is always
        // at inf which is undefined numerically
				conf .set_domain(dom).set(new_pos) = 1.;
				lapse.set_domain(dom).set(new_pos) = 1.;
      }
      // loop over all colocation points
		} while(new_pos.inc());
	} //end importing fields

  // Want to make sure there is no matter around the BH
  for(int d = space.BH; d < space.OUTER; ++d){
    logh.set_domain(d).annule_hard();
    phi.set_domain(d).annule_hard();
  }
  for(int d = space.BH; d < space.BH+2; ++d) {
    conf.set_domain(d).annule_hard();
    lapse.set_domain(d).annule_hard();
  }

  // employ standard spectral expansion, compatible with the given paraties
  conf.std_base();
	lapse.std_base();
	logh.std_base();
  shift.std_base();
  phi.std_base();

  // save everything to a binary file
  bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
}