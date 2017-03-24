/*
    Copyright 2017 Philippe Grandclement

    This file is part of Kadath.

    Kadath is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Kadath is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kadath.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __BIN_BH_HPP_
#define __BIN_BH_HPP_

#include "space.hpp"
#include "spheric.hpp"
#include "homothetic.hpp"
#include "bispheric.hpp"
namespace Kadath  {

/**
 * Spacetime intended for binary black hole configurations (see constructor for details about the domains used)
 * The radii of the holes are variables.
 * \ingroup domain
 */

class Space_bin_bh : public Space {
  
  protected:
    double a_minus ; ///< X-absolute coordinate of the center of the first sphere.
    double a_plus ; ///< X-absolute coordinate of the center of the second sphere.
  
     public:  

/**
     * Standard constructor    
     * @param ttype [input] : the type of basis.
     * @param dist [input] : distance \f$ d \f$ between the centers of the two spheres.
     * @param r1 [input] : radius \f$ r_1 \f$ of the first sphere.
     * @param r2 [input] : radius \f$ r_2 \f$ of the second sphere.
     * @param rext [input] : radius \f$ R \f$ of the outer boundary of the bispherical part.
     * @param nr [input] : number of points in each dimension 
     *
     * The various domains are then :
     * \li One \c Domain_nucleus, of radius \f$ r_1/2 \f$, centered at \f$ x_1 \f$.
     * \li One \c Domain_shell_outer_homothetic of with \f$ r_1/2 < r < r_1 \f$, centered at \f$ x_1 \f$.
     * \li One \c Domain_shell_inner_homothetic of with \f$ r_1 < r < 2 r_1 \f$, centered at \f$ x_1 \f$.
     * \li One \c Domain_nucleus, of radius \f$ r_2/2 \f$, centered at \f$ x_2 \f$.
     * \li One \c Domain_shell_outer_homothetic of with \f$ r_2/2 < r < r_2 \f$, centered at \f$ x_2 \f$.
     * \li One \c Domain_shell_inner_homothetic of with \f$ r_2 < r < 2 r_2 \f$, centered at \f$ x_2 \f$.
     * \li One \c Domain_bispheric_chi_first near the first sphere.
     * \li One \c Domain_bispheric_rect near the first sphere.
     * \li One \c Domain_bispheric_eta_first inbetween the two spheres.
     * \li One \c Domain_bispheric_rect near the second sphere.
     * \li One \c Domain_bispheric_chi_first near the second sphere.
     * \li One \c Domain_compact centered on the origin, with inner radius \f$ R \f$.
     */
	Space_bin_bh(int ttype, double dist, double rbh1, double rbh2, double rext, int nr) ;
	
	Space_bin_bh (FILE*) ; ///< Constructor from a file
     
        virtual ~Space_bin_bh() ; ///< Destructor
	virtual void save(FILE*) const ;

	virtual int nbr_unknowns_from_variable_domains() const ;
	virtual void affecte_coef_to_variable_domains(int& , int, Array<int>&) const ;
	virtual void xx_to_ders_variable_domains(const Array<double>&, int&) const ;
	virtual void xx_to_vars_variable_domains(System_of_eqs*, const Array<double>&, int&) const ;
	virtual Array<int> get_indices_matching_non_std(int dom, int bound) const ;
	
	/**
	* Adds a bulk equation and two matching conditions.
  	* The two nucleii and inner shells are not part of the computational space.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the bulk equation.
	* @param rac : the string describing the first matching condition.
	* @param rac_der : the string describing the second matching condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_eq (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der, int nused=-1, Array<int>** pused=0x0)  ; 
	/**
	* Sets a boundary condition at the inner radius of the first outer shell.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the boundary condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_bc_sphere_one (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  ;
	/**
	* Sets a boundary condition at the inner radius of the second outer shell.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the boundary condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_bc_sphere_two (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  ;
	/**
	* Sets a boundary condition at the outer boundary.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the boundary condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_bc_outer (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  ;
	 
	/**
	* Adds an equation being a surface integral at infinity.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the equation (should contain something like integ(f)=b)
	*/
	void add_eq_int_inf (System_of_eqs& syst, const char* eq) ;
 	/**
	* Adds an equation being a surface integral on the first sphere.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the equation (should contain something like integ(f)=b)
	*/
	void add_eq_int_sphere_one (System_of_eqs& syst, const char* eq) ;
	/**
	* Adds an equation being a surface integral on the first sphere.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the equation (should contain something like integ(f)=b)
	*/
	void add_eq_int_sphere_two (System_of_eqs& syst, const char* eq) ;
	/**
	* Adds an equation saying that one coefficient of a field is zero (at infinity)
	* @param syst : the \c System_of_eqs.
	* @param f : the field
	* @param jtarget : the index \f$\theta\f$ of the mode that must vanished.
	* @param ktarget : the index \f$\varphi\f$ of the mode that must vanished.
	*/
	void add_eq_zero_mode_inf (System_of_eqs& syst, const char* f, int jtarget, int ktarget) ;
} ;
}
#endif
