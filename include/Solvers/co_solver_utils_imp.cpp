#include "EOS/EOS.hh"
#include "coord_fields.hpp"
#include "ns_3d_xcts/ns_3d_xcts_solver.hpp"
#include "bco_utilities.hpp"

/**
 * \addtogroup Solver_utils
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
  
template <std::size_t s_type,typename config_t>
void setup_co(config_t& bconfig) {
  auto& fields = bconfig.return_fields();

  int type_coloc = CHEB_TYPE;
  Dim_array res(bconfig(DIM));
  res.set(0) = bconfig(BCO_PARAMS::BCO_RES);
  res.set(1) = bconfig(BCO_PARAMS::BCO_RES);
  res.set(2) = bconfig(BCO_PARAMS::BCO_RES)-1;

  Point center(bconfig(BCO_PARAMS::DIM));
  for (int i = 1; i <= bconfig(BCO_PARAMS::DIM); i++)
    center.set(i) = 0;
  
  const int shells = (int)bconfig(BCO_PARAMS::NSHELLS);
  int ndom = 4 + shells;
  std::vector<double> bounds(ndom - 1);
  if constexpr (s_type == NODES::BH){
    // Estimate Radius of BH based on Schwarzschild radius and an 
    // estimate for Psi on the horizon based on prev. BH solutions
    if(!bconfig.control(CONTROLS::USE_CONFIG_VARS)) {
      bconfig.set(BCO_PARAMS::RIN) = bconfig(BCO_PARAMS::MCH) * Kadath::bco_utils::invpsisq;
      bconfig.set(BCO_PARAMS::RMID) = 2 * bconfig(BCO_PARAMS::MCH) * Kadath::bco_utils::invpsisq;
      bconfig.set(BCO_PARAMS::ROUT) = 4 * bconfig(BCO_PARAMS::RMID);
    }    
    Kadath::bco_utils::set_isolated_BH_bounds(bounds, bconfig);
    #ifdef DEBUG
    Kadath::bco_utils::print_bounds("BH", bounds);
    #endif
    Space_adapted_bh space(type_coloc, center, res, bounds);

    write_bh_init_setup_tofile_XCTS(space, bconfig);
  }
  else if constexpr (s_type == NS) {
    const double h_cut = bconfig.template eos<double>(EOS_PARAMS::HCUT);
    const std::string eos_file = bconfig.template eos<std::string>(EOS_PARAMS::EOSFILE);
    const std::string eos_type = bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE);
    
    // Lambda to reduce code based on EOSTYPE
    auto gen_NS = [&](auto tov) {
      Kadath::bco_utils::set_NS_bounds(bounds, bconfig);
      
      // generate a full single star space including compactification to infinity
      Space_spheric_adapted space(type_coloc, center, res, bounds);
      
      write_ns_init_setup_tofile_XCTS(space, bconfig, *tov);
    };
    
    if(!bconfig.control(CONTROLS::USE_CONFIG_VARS)) {
      if(eos_type == "Cold_PWPoly") {
        using eos_t = ::Kadath::Margherita::Cold_PWPoly;
        EOS<eos_t, eos_var_t::PRESSURE>::init(eos_file, h_cut);
        auto tov = setup_ns_config_from_TOV<eos_t>(bconfig);
        gen_NS(std::move(tov));
      } else if(eos_type == "Cold_Table") {
        using eos_t = ::Kadath::Margherita::Cold_Table;

        const int interp_pts = (bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS) == 0) ? \
                                2000 : bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS);

        EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
        auto tov = setup_ns_config_from_TOV<eos_t>(bconfig);
        gen_NS(std::move(tov));
      }
      else { 
        std::cerr << eos_type << " is not recognized.\n";
        std::_Exit(EXIT_FAILURE);
      }
    }
  }
  else 
    std::cerr << "BCO type not implemented. \n";
  
}

template<typename config_t>
void write_bh_init_setup_tofile_XCTS(Space_adapted_bh& space, config_t& bconfig) {
  Base_tensor basis(space, CARTESIAN_BASIS);

  // setup fields
  Scalar lapse(space);
  lapse = 1.;
  lapse.set_domain(0).annule_hard();
  lapse.set_domain(1).annule_hard();

  Scalar conf(lapse);
  
  // set a better estimate for PSI on the horizon
  conf.set_domain(2) = Kadath::bco_utils::psi;
  conf.std_base();
  lapse.std_base();

  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();
  shift.std_base();
  // end setup fields
  
  Kadath::bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
}

template<typename tov_t, typename config_t>
void write_ns_init_setup_tofile_XCTS(Space_spheric_adapted& space, config_t& bconfig, tov_t& tov) {
  using eos_t = typename tov_t::eos_t;
  enum ltpQ { LAPSE=0, RHO, CONF }; 
  const int ndom = space.get_nbr_domains();
  Base_tensor basis(space, CARTESIAN_BASIS);
  // setup fields
  Scalar lapse(space);
  lapse = 1.;

  Scalar conf(lapse);
  
  // H = log(h), the logarithm of the specific enthalpy
  Scalar logh(space);
  logh.annule_hard();
  
  // interpolate TOV solution
  auto lintp = setup_interpolator_from_TOV(tov); 
  
	// generate a radius field for use with the 1D interpolators
  CoordFields<Space_spheric_adapted> cfg(space);
  Scalar r_field(cfg.radius());

  // update the fields based on TOV solution for a given domain
  auto update_fields= [&](const size_t dom) {
    Index pos(space.get_domain(dom)->get_nbr_points());
    do {
      double rval = r_field(dom)(pos);
      auto all_ltp = lintp.interpolate_all(rval);
      auto rho = (all_ltp[ltpQ::RHO] <= 0) ? 1e-15 : all_ltp[ltpQ::RHO];

      auto h = EOS<eos_t, eos_var_t::DENSITY>::h_cold__rho(rho);
      if(pos(0) == 0 && pos(1) == 0 && pos(2) == 0)
        bconfig.set(HC) = h;
      logh.set_domain(dom).set(pos) = (h < 1) ? 0. : std::log(h); 
      lapse.set_domain(dom).set(pos) = all_ltp[ltpQ::LAPSE];
      conf.set_domain(dom).set(pos) = all_ltp[ltpQ::CONF];
    }while(pos.inc());
  };
  for(int i = 0; i < 3; ++i)
    update_fields(i);
  
  for (int d = 2; d < ndom; ++d)
    logh.set_domain(d).annule_hard();

  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();
  shift.std_base();
  logh.std_base();
  conf.std_base();
  lapse.std_base();
  // end setup fields
  
  Kadath::bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh);
}

template<typename eos_t, typename config_t>
auto setup_ns_config_from_TOV(config_t& bconfig) {
  using namespace Kadath::Margherita;
  auto tov = std::make_unique<MargheritaTOV<eos_t>>();
  bool use_Mmax = tov->solve_for_MADM(bconfig(BCO_PARAMS::MADM));
  
  if(use_Mmax)
    bconfig.set(BCO_PARAMS::MADM) = tov->mass;
  
  bconfig.set(BCO_PARAMS::NC) = tov->rhoc;
  bconfig.set(BCO_PARAMS::HC) = EOS<eos_t, eos_var_t::PRESSURE>::h_cold__rho(bconfig(BCO_PARAMS::NC));
  bconfig.set(BCO_PARAMS::MB) = tov->baryon_mass;

  // update surface radius estimate
  bconfig.set(BCO_PARAMS::RMID) = tov->radius;
  bconfig.set(BCO_PARAMS::RIN) = 0.5 * bconfig(BCO_PARAMS::RMID);
  bconfig.set(BCO_PARAMS::ROUT) = 2.0 * bconfig(BCO_PARAMS::RMID);

  return std::move(tov);
}

template<typename tov_t>
auto setup_interpolator_from_TOV(tov_t& tov) {
  const size_t max_iter = tov.state.size();
  
  std::unique_ptr<double[]> radius_lin_ptr{new double[max_iter]};
  std::unique_ptr<double[]> conf_lin_ptr{new double[max_iter]};
  std::unique_ptr<double[]> lapse_lin_ptr{new double[max_iter]};
  std::unique_ptr<double[]> rho_lin_ptr{new double[max_iter]};
  
  for(auto j=0; j < max_iter; ++j) {
    //lapse est.
    auto phi = tov.state[j][tov.PHI];
    
    lapse_lin_ptr[j] = std::exp(phi);
    
    conf_lin_ptr[j] = tov.state[j][tov.CONF];

    radius_lin_ptr[j] = tov.state[j][tov.RISO];

    rho_lin_ptr[j] = tov.state[j][tov.RHOB];
  }
  
  // linear interpolation of conf(r_isotropic), lapse(r_isotropic), and rho(r_isotropic)
  // The order here dictates the enum ltpQ enum
  linear_interp_t<double,3> ltp(max_iter, std::move(radius_lin_ptr), std::move(lapse_lin_ptr), std::move(rho_lin_ptr), std::move(conf_lin_ptr)); 
  return ltp; 
}
}}
/** @}*/