/*
 * =====================================================================================
 *
 *       Filename:  tov.hh
 *
 *    Description:  Margherita TOV: Simple TOV solving routines
 *
 *        Version:  2.0
 *        Created:  1/07/2021
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Samuel David Tootle, tootle@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *   Notes: this is a minimal rewrite of the original TOV solver by ERM
 *
 * =====================================================================================
 */
#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <functional>

namespace Kadath {
namespace Margherita {
template <typename EOS> class MargheritaTOV {
  public:
  enum {RADIUS=0, RHOB, PRESS, MASSR, MASSB, PHI, RISO, CONF, NUMVAR};
  typedef EOS eos_t;

  private:
  using ary_t = std::array<double,NUMVAR>;
  using state_t = std::vector<ary_t>;
  
  static constexpr size_t max_iter = 10000;
  static constexpr double dr0 = 10. / 10000;
  const double loop_eps = 1e-6;

  public:
  bool adaptive{true};
  bool rk45{true};
  double radius{0};
  double arealr{0};
  double mass{0};
  double rhoc{0};
  double baryon_mass{0};
  state_t state{};

  void reset(const double& rho_init) {
    state.clear();
    state.reserve(max_iter);
    radius = 0;
    arealr = 0;
    mass = 0;
    baryon_mass = 0;
    rhoc = rho_init;
  }
  
  /**
   * gen_initial
   *
   * here the initial conditions are set.  Update this function
   * in the event more variables are added.
   */
  inline ary_t gen_initial() {
    using namespace Kadath::Margherita;
    double eps;
    typename EOS::error_t err;
    double press = EOS::press_cold_eps_cold__rho(eps, rhoc, err);
    ary_t initial_conditions {};

    initial_conditions[RHOB] = rhoc;
    initial_conditions[PRESS] = press;
    initial_conditions[PHI] = 1.;
    return std::move(initial_conditions);
  }

  /** 
   * evolve
   *
   * update function for dy's - returns k * h
   * if more variables are added, one needs to add the relevant
   * equations here and update the class enum.
   *
   * @param f your input array of dy's
   * @param h step size
   */
  inline ary_t evolve(const ary_t& f, const double& h) const {
    typename EOS::error_t err;
    if(f[RADIUS] == 0) return f;

    ary_t res{};
    res.fill(0);
    auto m = f[MASSR];
    double r = f[RADIUS];
    double r2 = r * r;
    double r3 = r2* r;
    double dedp, press, rho, rhoE;
    press = f[PRESS];
    rho = EOS::rho_energy_dedp__press_cold(rhoE, dedp, press, err);
    
    // eq(178) Tichy. https://arxiv.org/abs/1610.03805
    res[MASSR] = 4. * M_PI * r2 * rhoE * h;
    // eq(180) Tichy. https://arxiv.org/abs/1610.03805
    res[PHI] = (m + 4. * M_PI * r3 * press) / (r * (r - 2 * m)) * h;
    // eq(179) Tichy. https://arxiv.org/abs/1610.03805
    res[PRESS] = -res[PHI] * (rhoE + press);

    double sch_fac = r / (r - 2. * m) ;
    // eq(18) integral only
    res[RISO] = (std::sqrt(sch_fac) - 1.) / r * h;
    res[MASSB] = std::sqrt(sch_fac) * res[MASSR];
    return res;
  }
  
  // RK45 with option for adaptive step sizes - enabled by default
  inline ary_t rk45_step(const ary_t& input, double& dr) const { 
    using namespace Kadath::Margherita;
    typename EOS::error_t err;
    auto eps = 5.e-8 * input[PRESS];

    auto k1 = evolve(input, dr);

    auto s = input;
    
    s[RADIUS] += (2./9.) * dr ;
    for(auto i = 1; i < NUMVAR; ++i) {
      s[i] += (2./9.) * k1[i];
    }
    
    auto k2 = evolve(s, dr);
    
    s = input;
    s[RADIUS] += dr / 3.;
    for(auto i = 1; i < NUMVAR; ++i) {
      s[i] += 1./12. * k1[i] + 0.25 * k2[i] ;
    }

    auto k3 = evolve(s, dr);
    
    s = input;
    s[RADIUS] += 0.75 * dr ;
    for(auto i = 1; i < NUMVAR; ++i) {
      s[i] += (69./128.) * k1[i] + (-243./128.) * k2[i] + (135./64.) * k3[i];
    } 

    auto k4 = evolve(s, dr);
    
    s = input;
    s[RADIUS] += dr ;
    for(auto i = 1; i < NUMVAR; ++i) {
      s[i] += (-17./12.) * k1[i] + (27./4.) * k2[i] + (-27./5.) * k3[i] + (16./15.) * k4[i];
    } 
    auto k5 = evolve(s, dr);

    s = input;
    s[RADIUS] += (5./6.) * dr ;
    for(auto i = 1; i < NUMVAR; ++i) {
      s[i] += (65./432.) * k1[i] + (-5./16.) * k2[i] + (13./16.) * k3[i] + (4./27.) * k4[i] + (5./144.) * k5[i];
    } 
    auto k6 = evolve(s, dr);
    
    auto res = input;
    for(int i = 1; i < NUMVAR; ++i)
      res[i] += (1./450.) * (47. * k1[i] + 216. * k3[i] + 64. * k4[i] + 15. * k5[i] + 108. * k6[i]);
    
    // manual update of tracked variables
    res[RADIUS] += dr;
    res[RHOB] = EOS::rho__press_cold(res[PRESS],err);
    
    ary_t error;
    for (int nn = 1; nn < NUMVAR; ++nn) {
      error[nn] = 1. / 300. *
                  (-2. * k1[nn] + 9. * k3[nn] - 64. * k4[nn] - 15. * k5[nn] +
                   72. * k6[nn]);
    };

    auto adaptive_step = [&]() {
      // Is the error too large?
      if (std::abs(error[PRESS]) > eps && state.size() > 2) {
        if (dr > 2. * dr0) {
          dr = dr * 0.5;
        } else {
          dr *= 0.9 * std::pow(eps / std::abs(error[PRESS]), 1./5.);
        }
        res = rk45_step(input, dr);
      } else if (std::abs(error[PRESS]) < 1. / 20. * eps) {
        dr *= 2.;
      }
    };
    if(adaptive)
      adaptive_step();
       
    return std::move(res);
  }
  
  // basic RK4 - mainly for testing
  inline ary_t rk_step(const ary_t& input, double& dr) const { 
    using namespace Kadath::Margherita;
    typename EOS::error_t err;
    auto eps = 5.e-8 * input[PRESS];

    auto k1 = evolve(input, dr);

    auto s1 = input;
    
    s1[RADIUS] += dr / 2.;
    for(auto i = 1; i < NUMVAR; ++i)
      s1[i] += k1[i] / 2.;
    
    auto k2 = evolve(s1, dr);
    
    auto s2 = input;
    s2[RADIUS] += dr / 2.;
    for(auto i = 1; i < NUMVAR; ++i)
      s2[i] += k2[i] / 2.;

    auto k3 = evolve(s2, dr);
    
    auto s3 = input;
    s3[RADIUS] += dr;
    for(auto i = 1; i < NUMVAR; ++i)
      s3[i] += k3[i];
    
    auto k4 = evolve(s3, dr);
    
    auto res = input;
    for(int i = 1; i < NUMVAR; ++i)
      res[i] += 1./6. * (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]);
    
    // manual update of tracked variables
    res[RADIUS] += dr;
    res[RHOB] = EOS::rho__press_cold(res[PRESS],err);
   
    return std::move(res);
  }

  // correct PHI based on enforcing Schwarzschild BC
  void correct_phi() {
    // eq. (5) Baumgarte Notes
    auto get_correction = [&](double m, double r, double phi) {
      return 0.5 * std::log(1. - 2 * m / r) - phi;
    };
    
    const double correction = 
      get_correction(state.back()[MASSR], state.back()[RADIUS], state.back()[PHI]);

    for(auto& el : state) 
      el[PHI] += correction;
  }

  /**
   * correct_isotropic_r
   * calculate the final isotropic radius and correct it based
   * on enforcing the Schwarzschild BC
   */
  inline void correct_isotropic_r() {
    // calculate isotropic radius
    // eq(18) from Baumgarte notes
    for(auto& el : state)  
      el[RISO] = el[RADIUS] * std::exp(el[RISO]);
      
    // eq.(19) from Baumgarte notes
    auto get_factor = [&](double m, double r, double rhat) {
      double f = std::sqrt(r*r - 2. * m * r) + r - m;
      f *= 0.5;
      f /= rhat;
      return f;
    };
    
    // get correction factor
    const double factor = 
      get_factor(state.back()[MASSR], state.back()[RADIUS], state.back()[RISO]);
    
    // apply correction
    for(auto& el : state) 
      el[RISO] *= factor;
  }

  /**
   * calculate_conformal_factor
   *
   * calculate conformal factor CONF such that
   * \gamma_ij = CONF^4 \eta_ij
   * eq(22) from Baumgarte notes
   */
  inline void calculate_conformal_factor() {
    for(auto& el : state) {
      const auto& r = el[RADIUS];
      const auto& rhat = el[RISO];
      el[CONF] = std::exp(0.5 * std::log(r / rhat));
    }
  }

  /**
   * solve
   *
   * find a solution to the TOV equations for the given central and EOS
   *
   * @param rho_init central density input
   * @param eps precision desired for finding ADM mass
   */
  inline void solve(const double rho_init, const bool fin = true) { 
    using std::placeholders::_1;
    using std::placeholders::_2;
    // determines which RK function to use instead of doing the
    // logic check multiple times in the loop.
    std::function<ary_t(ary_t&, double&)> rk_iteration = (rk45) ? 
      std::bind(&MargheritaTOV::rk45_step, this, _1, _2) : 
      std::bind(&MargheritaTOV::rk_step, this, _1, _2);
    
    // reset all member variables to defaults
    reset(rho_init);

    // initialize our starting dr in case adaptive RK is used
    double dr = dr0;

    // add initial conditions
    state.push_back(gen_initial());
    
    size_t count = 0;
    while(true && count < max_iter) {
      auto res = rk_iteration(state.back(), dr);
      
      double delm = std::abs(res[MASSR] - state.back()[MASSR]);
      if(delm < loop_eps && state.back()[PRESS] < loop_eps)
        break;
      state.push_back(res);
      count++;
    }
    if(count >= max_iter) {
      std::cerr << "Maximum iterations (" << max_iter << ") reached.\n"
        "The central density provided may not be able to produce a stable\n"
        "static TOV solution\n\n";
      std::_Exit(EXIT_FAILURE);
    }

    // erase initial conditions if a stable solution was found at all
    // if this is neglected, it causes issues in calculating RISO since
    // you get division by 0
    if(count > 0)
      state.erase(state.begin());
    else
      return;
   
    // save some computation if we just need ArealR, M, and Mb
    if(fin) {
      correct_phi();
      correct_isotropic_r();
      calculate_conformal_factor();
      radius = state.back()[RISO];
    }
    arealr = state.back()[RADIUS];
    mass = state.back()[MASSR];
    baryon_mass = state.back()[MASSB];
  }
  
  /**
   * solve_for_MADM
   *
   * Root finder that takes an Madm value and desired precision and attempts
   * to find a solution to the TOV equations for the given Madm and EOS by
   * varying the central density
   *
   * @param M_fin desired final ADM mass
   * @param eps precision desired for finding ADM mass
   */
  inline void solve_for_MADM(const double M_fin, const double eps = 1e-3) {
    double rhoL;
    const double rho_max0 = 1e-2;
    const double rho_min0 = 1e-4;
    double rho_min = rho_min0;
    double rho_max = rho_max0;
    double delrho = (rho_max - rho_min) / 100.;
  
    // find bracketing range based on the central density to 
    // to make sure a root exists - i.e. Madm can be found
    for(rhoL = rho_min0; rhoL < rho_max0; rhoL += delrho) {
      solve(rhoL, false);
      if(mass < M_fin)
        rho_min = rhoL;
      else if(mass > M_fin) {
        rho_max = rhoL;
        break;
      }
    }
    
    if(rhoL >= rho_max0) {
      std::cerr << "Maximum density (" << rhoL << ") reached.\n"
        "The density bracketing range may not be sufficiently constrained for the given EOS\n";
      std::_Exit(EXIT_FAILURE);
    }

    // from the discovered bracketing range, we now find the ADM Mass using
    // a basic bisection of the central density
    size_t count = 0;
    double diff = 1;
    while( std::fabs(diff) > eps && count < max_iter) {
      rhoL = (rho_max + rho_min) / 2.;      
      solve(rhoL, false);
      diff = (mass-M_fin);
      
      #ifdef margerita_check_all
      std::cout << "M: " << mass
                << "\nMb: " << baryon_mass
                << "\nArealR: " << arealr
                << "\nNc: " << rhoL
                << "\nDiff: " << diff 
                << "\nIter: " << count << "\n\n";
      #endif
      if( diff > 0) 
        rho_max = rhoL;
      else
        rho_min = rhoL;
      count++;
    }
    if(count >= max_iter) {
      std::cerr << "Maximum iterations (" << max_iter << ") reached.\n"
        "The mass may be too high or the range of rho may not be sufficiently wide\n";
      std::_Exit(EXIT_FAILURE);
    }
      
    correct_phi();
    correct_isotropic_r();
    calculate_conformal_factor();
    radius = state.back()[RISO];
  }
};
}}
