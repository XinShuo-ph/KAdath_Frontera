/*
    Copyright 2022 Philippe Grandclement

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

#ifndef __ADAPTED_HOMOTHETIC_HPP_
#define __ADAPTED_HOMOTHETIC_HPP_

#include "space.hpp"
#include "term_eq.hpp"
#include "homothetic.hpp"
#include "metric.hpp"

namespace Kadath {

/**
 * The \c Space_spheric_adapted _homothetic class fills the space with one nucleus, one spherical shell of varying outer radius, one spherical shell of varying inside radius , several standard shells and a compactified domain, all centered on the same point.
 * \ingroup domain
 */
class Space_spheric_adapted_homothetic : public Space {
     public:
	/**
     	* Standard constructor ; all the shells are initially not deformed.
	* @param ttyp [input] : the type of basis.
	* @param cr [input] : absolute coordinates of the center.
	* @param nbr [input] : number of points in each domain.
	* @param bounds [input] : radii of the various shells (and also determines the total number of domains).
	*/
	Space_spheric_adapted_homothetic  (int ttyp, const Point& cr, const Dim_array& nbr, const Array<double>& bounds) ;
	
	Space_spheric_adapted_homothetic (FILE*) ; ///< Constructor from a file
	virtual ~Space_spheric_adapted_homothetic () ; ///< Destructor        
	virtual void save(FILE*) const ;
	
	virtual int nbr_unknowns_from_variable_domains() const ;
	virtual void affecte_coef_to_variable_domains(int& , int, Array<int>&) const ;
	virtual void xx_to_ders_variable_domains(const Array<double>&, int&) const ;
	virtual void xx_to_vars_variable_domains(System_of_eqs*, const Array<double>&, int&) const ;

	virtual Array<int> get_indices_matching_non_std(int, int) const ;
} ;
}
#endif
