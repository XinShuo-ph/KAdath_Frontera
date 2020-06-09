/*
    Copyright 2020 Philippe Grandclement

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

#ifndef __BBH_HPP_
#define __BBH_HPP_

#include "adapted.hpp"
#include "bispheric.hpp"
namespace Kadath  {

/**
 * Spacetime intended for binary black hole configurations in full general relativity (see constructor for details about the domains used)
 * The shape of the horizons are variables.
 * \ingroup domain
 */

class Space_bbh : public Space {
  
     public:  

/**
     * Standard constructor    
     * @param ttype [input] : the type of basis.
     * @param dist [input] : distance between the centers of the two spheres.
     * @param r1 [input] : radius of the first hole (assumed spherical at this point).
     * @param r2 [input] : radius of the second hole (assumed spherical at this point).
     * @param rbi [input] : radius of the outer boundary of the bispherical part.
     * @param rext [input] : outer radius of the. 
     * @param nr [input] : number of points in each dimension 
     *
     * The various domains are then :
     * \li One \c Domain_shell_outer_adapted around the first hole.
     * \li One \c Domain_shell_inner_adapted around the first hole.    
     * \li One \c Domain_shell_outer_adapted around the second hole.
     * \li One \c Domain_shell_inner_adapted around the second hole.
     * \li One \c Domain_bispheric_chi_first near the first hole.
     * \li One \c Domain_bispheric_rect near the first hole.
     * \li One \c Domain_bispheric_eta_first inbetween the two holes.
     * \li One \c Domain_bispheric_rect near the second hole.
     * \li One \c Domain_bispheric_chi_first near the second hole.
     * \li One \c Domain_shell centered on the origin, ecnompassing the two holes.
     */
	Space_bbh(int ttype, double dist, double rbh1, double rbh2, double rbi, double rext, int nr) ;
	
	Space_bbh (FILE*) ; ///< Constructor from a file
     
	/**
	* Sets a matching condition accross all the bispheric domain (intended for a second order equation).
	* @param syst : the \c System_of_eqs.
	* @param rac : the string describing the boundary condition.
	* @param list : list of the components to be considered.
	*/
	void add_matching (System_of_eqs& syst, const char* rac, const List_comp& list)  ;

	/**
	* Sets a matching condition accross all the bispheric domain (intended for a second order equation).
	* @param syst : the \c System_of_eqs.
	* @param rac : the string describing the boundary condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_matching (System_of_eqs& syst, const char* rac, int nused=-1, Array<int>** pused=0x0)  ;
	

        virtual ~Space_bbh() ; ///< Destructor
	virtual void save(FILE*) const ;

	virtual int nbr_unknowns_from_variable_domains() const ;
	virtual void affecte_coef_to_variable_domains(int& , int, Array<int>&) const ;
	virtual void xx_to_ders_variable_domains(const Array<double>&, int&) const ;
	virtual void xx_to_vars_variable_domains(System_of_eqs*, const Array<double>&, int&) const ;
	virtual Array<int> get_indices_matching_non_std(int dom, int bound) const ;
} ;
}
#endif
