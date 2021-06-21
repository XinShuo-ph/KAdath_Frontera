/*
 * This file is part of the KADATH library.
 * Author: Samuel Tootle
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
#include <sstream>
using namespace Kadath;

#define EQUI -11

/**
 * @namespace bco_utils
 * Kadath specific utilities to reduce duplicate code in 
 * initial data solvers for compact objects and make code (hopefully)
 * easier to read.
 */
namespace bco_utils {
// approximate value of psi on the horizon
// to make coordinate estimates
constexpr double psi = 1.55 ;
constexpr double psisq = psi * psi;
constexpr double invpsisq = 1. / psisq;

/**
 * bco_utils::KadathPNOrbitalParams \n 
 * 3.5 PN orbital parameters for irrotational compact object
 * This is a very rough, but useful guess for any CO
 *
 * @param[input]  bconfig:  binary Configurator file
 * @param[input]  Madm1: ADM mass (FIXED_MADM or MCH)
 * @param[input]  Madm2: ADM mass (FIXED_MADM or MCH)
 */
void KadathPNOrbitalParams(kadath_config_boost<BIN_INFO>& bconfig, double Madm1, double Madm2);

/**
 * bco_utils::save_to_file\n 
 * Standard print config, space, and fields to their respective files
 * Config is saved to base_fname.info, space and fields in base_fname.dat
 * Output directory is stored in the config file. See config_utils_boost.hpp
 *
 * @tparam config_t Configurator type
 * @tparam space_t Space type
 * @tparam fields_t parameter pack of Kadath field types to save to file
 * @param[input]  base_fname: base filename (no extension)
 * @param[input]  bconfig:  Configurator file
 * @param[input]  space:  Numerical space
 * @param[input]  fields: Parameter pack of fields to save to file
 */
template <typename space_t, typename config_t, typename... fields_t>
void save_to_file(std::stringstream& base_fname, space_t& space, config_t& bconfig, fields_t&... fields){
  bconfig.set_filename(base_fname.str());
  bconfig.write_config();
  std::string kadath_filename = bconfig.space_filename();
  FILE *ff = fopen(kadath_filename.c_str(), "w");
  space.save(ff) ;

  //All we want to do is call fields.save(ff), but C++ is picky on how to do that with packs
  (void) std::initializer_list<int>{(fields.save(ff),0)...};
  fclose(ff);
}

/**
 * bco_utils::save_to_file
 * 
 * Specialization that excludes the base_filename stringstream.\n 
 * Standard print config, space, and fields to their respective files.\n 
 *  - Config is saved to its set preset filename\n 
 *  - Fields are save to the paired filename. See config_utils_boost.hpp.\n 
 *
 * @tparam config_t Configurator type
 * @tparam space_t Space type
 * @tparam fields_t parameter pack of Kadath field types to save to file
 * @param[input]  space:  Numerical space
 * @param[input]  bconfig:  Configurator file
 * @param[input]  fields: Parameter pack of fields to save to file
 */
template <typename space_t, typename config_t, typename... fields_t>
void save_to_file(space_t& space, config_t& bconfig, fields_t&... fields) {
  bconfig.write_config();
  
  std::string kadath_filename = bconfig.space_filename();
  FILE *ff = fopen(kadath_filename.c_str(), "w");
  space.save(ff) ;

  //All we want to do is call fields.save(ff), but C++ is picky on how to do that with packs
  (void) std::initializer_list<int>{(fields.save(ff),0)...};
  fclose(ff);
}

/**
 * bco_utils::get_radius
 * 
 * When all we want is a simply radius value (inner or outer radius)
 *
 * Note: For variable domains, this is just R at a single point on the surface
 *
 * @tparam dom_t Domain type
 * @param[input]  dom: pointer to a domain with type dom_t
 * @param[input]  bound: need to know if you want inner or outer radius
 * @return a radius of the domain.   */
template<typename dom_t>
double get_radius (const dom_t* dom, const int bound) {
    Index pos(dom->get_nbr_points());
    switch(bound) {
      case INNER_BC:
        break;
      case OUTER_BC:
        pos.set(0) = dom->get_nbr_points()(0) - 1;
        pos.set(1) = dom->get_nbr_points()(1) - 1;
        pos.set(2) = dom->get_nbr_points()(2) - 1;
        break;
      case EQUI:
        pos.set(0) = dom->get_nbr_points()(0) - 1;
        break;
      default:
        std::cout << "Unknown bound sent to get_radius: " << bound << std::endl;
        break;
    }
    return dom->get_radius()(pos);
}

/**
 * bco_utils::get_min_max
 * 
 * For variable domains, find the min and max of a Scalar field on a boundary
 *
 * @tparam space_t type of numerical space
 * @param[input]  field: Scalar field
 * @param[input]  dom: domain to find min and max of
 * @return {fmin, fmax} std::array containing fmin and fmax
 */
inline std::array<double,2> get_field_min_max(const Scalar& field, const int dom, const int bound = OUTER_BC) {
  const int npts_r = field.get_domain(dom)->get_nbr_points()(0);
  
  // position index to loop over
  Index pos(field.get_domain(dom)->get_nbr_points());
  // position index to update with a fixed radial boundary
  Index bpos(field.get_domain(dom)->get_nbr_points());

  int r_bound = 0;
  switch(bound) {
    case INNER_BC:
      break;
    case OUTER_BC:
      r_bound = npts_r - 1;
      break;
    default:
      std::cout << "Unknown bound sent to get_field_min_max: " << bound << std::endl;
      break;
  } 

	double fmax = field(dom)(pos);
  double fmin = field(dom)(pos);
  
  // FIXME this currently loops over theta and phi multiple times, but the syntax is straightforward.
  do {
      bpos.set(0) = r_bound;
      bpos.set(1) = pos(1);
      bpos.set(2) = pos(2);
  		double f = field(dom)(bpos);

  		if(f > fmax)
  			fmax = f;
  		if(f < fmin)
  			fmin = f;
  } while(pos.inc());
  return std::array<double,2>{fmin,fmax};
}

/**
 * bco_utils::get_rmin_rmax
 * 
 * For variable domains, find the min and max radii
 *
 * @tparam space_t type of numerical space
 * @param[input]  space: Numerical space
 * @param[input]  dom: domain to find min and max radii of
 * @return {rmin, rmax} std::array containing rmin and rmax
 */
template<typename space_t>
std::array<double,2> get_rmin_rmax(const space_t& space, const int dom) {
  const int npts_r = space.get_domain(dom)->get_nbr_points()(0);
  Index pos(space.get_domain(dom)->get_nbr_points());
  pos.set_start();
  pos.set(0) = npts_r - 1;

	double rmax = space.get_domain(dom)->get_radius()(pos);
  double rmin = space.get_domain(dom)->get_radius()(pos);
  do {
      pos.set(0) = npts_r - 1;
  		double r = space.get_domain(dom)->get_radius()(pos);

  		if(r > rmax)
  			rmax = r;
  		if(r < rmin)
  			rmin = r;
  } while(pos.inc());
  return std::array<double,2>{rmin,rmax};
}

/**
 * bco_utils::update_adapted_field
 *
 * extrapolate variable fields to domains that are/should be zero based on the coefficients on the boundary, e.g:
 * 	- extrapolate space-time fields into the horizon of the BH from outside the horizon
 * 	- extrapolate matter fields from inside the star to outside of the star
 * 
 * Note 1:\n 
 * This is important before importing fields from one numerical space to another as interpolating the fields to the new\n 
 * space may have problems if the adjacent domain is zero.  After import, remember to set relevant domains back to zero\n 
 * 
 * Note 2:\n 
 * could be generalized to Tensors.  Currently only done for Scalars
 *
 * @tparam dom_t type of domain.
 * @param[input]  res: reference to modifiable Scalar field that contains the current solution to the field.  This will be changed!
 * @param[input]  from_dom: index of domain to extrapolate FROM 
 * @param[input]  to_dom: index of domain to extrapolate TO
 * @param[input]  dom: pointer to domain to extrapolate TO
 * @param[input]  bound: which boundary we extrapolating from (OUTER_BC, INNER_BC)
 * @param[output] res: Modified scalar field with values interpolated from from_dom to to_dom
 */
template<typename dom_t>
void update_adapted_field (Scalar& res, const int from_dom, const int to_dom, const dom_t* dom, const int bound) {
  Tensor* ptr = &res;
  Array<int> doms(1);
  doms.set(0) = from_dom;
  Scalar import = dom->import(to_dom, bound, 1, doms, &ptr);
  ptr->set().set_domain(to_dom) = import(to_dom);
}

/**
 * bco_utils::interp_adapted_mapping\n 
 * Update the mapping of an adapted domain based on an old numerical space before importing fields
 *
 * @tparam adapted_t type of adapted domain
 * @param[input]  new_shell: const pointer to adapted domain of the new shell that we're creating a mapping for
 * @param[input]  old_outer_adapted_dom: index of old outer adapted dom that we are creating a mapping from
 * @param[input]  old_radius_field: Scalar field describing the radii of all the fields in the numerical space to interpolate on
 */
template<typename adapted_t>
void interp_adapted_mapping(const adapted_t* new_shell, const int old_outer_adapted_dom, const Scalar& old_radius_field) {
  Val_domain new_mapping = new_shell->get_radius();

  //Need to normalize by the constant radius at the INNER_BC of the outer_adapted shell in the old space
  auto old_shell = old_radius_field.get_space().get_domain(old_outer_adapted_dom);
  double rinner = get_radius(old_shell, INNER_BC);

  Index new_pos(new_shell->get_nbr_points());
  double xc_new = new_shell->get_center()(1);
  double xc_old = old_shell->get_center()(1);
  do {
    double x = new_shell->get_cart(1)(new_pos) - xc_new;
    double y = new_shell->get_cart(2)(new_pos);
    double z = new_shell->get_cart(3)(new_pos);

    double r = std::sqrt(x*x + y*y + z*z);

    x /= r / rinner;
    y /= r / rinner;
    z /= r / rinner;

    Point absol(3);
    absol.set(1) = x + xc_old;
    absol.set(2) = y;
    absol.set(3) = z;

    new_mapping.set(new_pos) = old_radius_field.val_point(absol);

  } while(new_pos.inc());

  new_mapping.std_base();
  new_shell->set_mapping(new_mapping);
}

/**
 * bco_utils::get_center\n 
 * Quickly get the X coordinate of the  cartesian center of a domain\n 
 * while using the correct Index and letting it be deleted when we're done.
 *
 * @tparam space_t type of numerical space
 * @param[input]  space: numerical space
 * @param[input]  dom: index of domain of interest
 * @return coordinate center on the X axis since we're always at (X,0,0)
 */
template<typename space_t>
double get_center (const space_t& space, const int dom) {
    Index pos(space.get_domain(dom)->get_nbr_points());
    return space.get_domain(dom)->get_cart(1)(pos);
}

/**
 * bco_utils::get_boundary_val\n 
 * Quickly get the value at a specified boundary of a scalar field in a given domain\n 
 * and letting the setup be deleted when we're done.
 *
 * @tparam space_t type of numerical space
 * @param[input]  space: numerical space
 * @param[input]  dom: index of domain of interest
 * @param[input]  field: scalar field of interest
 * @param[input]  bound: bondary of interest
 * @return value at the (0,0,0) collocation point
 */
inline double get_boundary_val (const int dom, const Scalar& field, const int bound=INNER_BC) {
  auto this_domain = field(dom).get_domain(); 
  Index pos(this_domain->get_nbr_points());
  switch(bound) {
    case INNER_BC:
      break;
    case OUTER_BC:
      pos.set(0) = this_domain->get_nbr_points()(0) - 1;
      pos.set(1) = this_domain->get_nbr_points()(1) - 1;
      pos.set(2) = this_domain->get_nbr_points()(2) - 1;
      break;
    case EQUI:
      pos.set(0) = this_domain->get_nbr_points()(0) - 1;
      break;
    default:
      std::cout << "Unknown bound sent to get_boundary_val: " << bound << std::endl;
      break;
  }
  return field(dom)(pos);
}

/**
 * bco_utils::set_radius\n 
 * Quickly get the basic outer radius of a given domain to update the corresponding\n 
 * radius in the config file based on the provides indexes 
 * 
 * @tparam space_t type of numerical space
 * @param[input]  dom: index of domain of interest
 * @param[input]  space: numerical space
 * @param[input]  bconfig: configuration file
 * @param[input]  idxs: 1 or more indexes related to the config file entry to update
 */
template<typename config_t, typename space_t, typename... idx_t>
void set_radius (const int& dom, const space_t& space, config_t& bconfig, const idx_t... idxs) {
    bconfig.set(idxs...) = get_radius(space.get_domain(dom), OUTER_BC);
}

/**
 * bco_utils::set_NS_bounds
 *
 * Set boundaries for a NS companion based on the config file parameters\n 
 * 
 * @tparam space_t type of numerical space
 * @param[input]  bounds: array of bounds to update
 * @param[input]  bconfig: configuration file
 * @param[input]  bco index: of the NS in the configuration file
 */
template<typename ary_t, typename config_t>
void set_NS_bounds (ary_t& bounds, config_t& bconfig, const int bco)  {
  //sizes
  const int size    = bounds.size();
  const int nshells = bconfig(NSHELLS, bco);

  //indexes
  const int rin     = 0;
  const int rout    = size - 1;
  const int r       = rin + nshells + 1;

  bounds[rin]       = bconfig(RIN, bco);
  bounds[r]         = bconfig(RMID, bco);
  bounds[rout]      = bconfig(ROUT, bco);

  double lower      = bounds[rin] + (bounds[r] - bounds[rin]) * 0.8;
  double delta_r    = (bounds[r] - lower) / (nshells + 1);

  for(int shell = rin+1; shell <= nshells; ++shell){
    bounds[shell] = delta_r * shell + lower;
  }
}

/**
 * bco_utils::set_isolated_BH_bounds
 *
 * Set boundaries for an isolated  BH based on the config file parameters
 * 
 * @tparam space_t type of numerical space
 * @param[input]  bounds: array of bounds to update
 * @param[input]  bconfig: configuration file
 */
template<typename ary_t, typename config_t>
void set_isolated_BH_bounds (ary_t& bounds, config_t& bconfig)  {

  // Configurator index of companion compact object
  const int size = bounds.size();
  
  //indexes
  const int rin  = 0;
  const int r    = rin + 1;
  int rout       = size - 1;

  bounds[rin]  = bconfig(RIN);
  bounds[r]    = bconfig(RMID);
  
  double delrBH = bconfig(ROUT) - bconfig(RMID);
  
  const int shells = bconfig(NSHELLS);
  double delr = delrBH / (bconfig(NSHELLS) + 1);

  for(int i = 0, b = r+1; i < shells; ++b, ++i)
    bounds[b] = bounds[r] + delr * (i + 1);
  
  bounds[rout] = bconfig(ROUT);
}

/**
 * bco_utils::set_BH_bounds
 *
 * Set boundaries for a BH companion based on the config file parameters
 * and it's companion
 * 
 * @tparam space_t type of numerical space
 * @param[input]  bounds: array of bounds to update
 * @param[input]  bconfig: configuration file
 * @param[input]  bco: index of the BH in the configuration file
 */
template<typename ary_t, typename config_t>
void set_BH_bounds (ary_t& bounds, config_t& bconfig, const int bco, const bool adapt_shells = false)  {

  // Configurator index of companion compact object
  const int bco2 = (bco == BCO1) ? BCO2 : BCO1;
  const int size = bounds.size();
  
  //indexes
  const int rin  = 0;
  const int r    = rin + 1;
  int rout       = size - 1;

  bounds[rin]  = bconfig(RIN, bco);
  bounds[r]    = bconfig(RMID, bco);
  
  // We compare resolution to the companion and add shells as needed
  // outside the horizon
  double delrC  = bconfig(ROUT, bco2)- bconfig(RMID, bco2);
  double delrBH = bconfig(ROUT, bco) - bconfig(RMID, bco);
  const int min_shells = int(std::floor(delrBH / delrC)) - 1;
  
  if(min_shells > bconfig(NSHELLS, bco) && adapt_shells) {
    std::cerr << "Number of shells for BH is too low.  Minimum required shells is: " << min_shells
              << ". Setting NSHELLS to the minimum and continuing. \n";
    rout += (min_shells - bconfig(NSHELLS,bco));
    bconfig.set(NSHELLS, bco) = min_shells;
    bounds.resize(size+min_shells);
  }
  
  //setup temporary BH config file to reuse isolated_BH code
  kadath_config_boost<BCO_BH_INFO> bhconfig;
  for(int i = 0; i < NUM_BCO_PARAMS; ++i) { bhconfig.set(i) = bconfig.set(i, bco); }
  set_isolated_BH_bounds(bounds, bhconfig);
}

/**
 * bco_utils::print_bounds
 *
 * helper function for printing domain boundaries in the readers - requires foreach compatible
 * container
 * 
 * @tparam space_t type of numerical space
 * @param[input]  bounds: array of bounds
 * @param[input]  name: string with the name you want associated with the bounds (e.g. NS1)
 */
template<typename ary_t>
void print_bounds(std::string name, const ary_t& bary) 
{
  std::cout << name << ": ";
  for(auto& e : bary)
  {
    std::cout << e << " ";
  }
  std::cout << std::endl;
 
}

/**
 * bco_utils::mirr_from_mch
 *
 * helper function to compute Mirr from MCH
 * 
 * @param[input]  chi: dimensionless spin parameter
 * @param[input]  mch: Christodoulou mass of the BH
 */
inline double mirr_from_mch(const double chi, const double mch){
  return std::sqrt((1+sqrt(1-chi*chi)) / 2.) * mch;
}

/**
 * bco_utils::syst_mch
 *
 * helper function to compute MCH from a system of equations
 * 
 * @param[input]  syst: system of equations
 * @param[input]  space: numerical space
 * @param[input]  eq: equation string to lookup for spin equation (e.g. Sint1, Sint2)
 * @param[input]  dom: domain index to evalute the integral (i.e. inner homothetic)
 */
template<typename space_t>
double syst_mch(System_of_eqs& syst, const space_t& space, const std::string eq, const int dom){
  Val_domain integS(syst.give_val_def(eq.c_str())()(dom));
  double S = space.get_domain(dom)->integ(integS, INNER_BC);
  
  Val_domain integMsq(syst.give_val_def("intMsq")()(dom));
  double Mirrsq = space.get_domain(dom)->integ(integMsq, INNER_BC);
  double Mirr = std::sqrt(Mirrsq);
  return std::sqrt(Mirrsq + S * S / 4. / Mirrsq );
}

}
