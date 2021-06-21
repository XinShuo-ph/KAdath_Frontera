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

#include "bispheric.hpp"
#include "val_domain.hpp"
#include "array.hpp"

namespace Kadath{
int mult_cos_1d (int, Array<double>&) ;
int mult_sin_1d (int, Array<double>&) ;
int div_x_1d (int, Array<double>&) ;

Val_domain Domain_bispheric_chi_first::mult_cos_phi (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base.allocate (nbr_coefs) ;	
	*res.base.bases_1d[2] = *so.base.bases_1d[2] ;
	
	Index index (res.base.bases_1d[1]->get_dimensions()) ;

	// Inversion in theta :
	do {
		switch ((*so.base.bases_1d[1])(index)) {
			case CHEB_EVEN :
				res.base.bases_1d[1]->set(index) = CHEB_ODD ;
				break ;
			case CHEB_ODD :
				res.base.bases_1d[1]->set(index) = CHEB_EVEN ;
				break ;
			case LEG_ODD :
				res.base.bases_1d[1]->set(index) = LEG_EVEN ;
				break ;
			case LEG_EVEN :
				res.base.bases_1d[1]->set(index) = LEG_ODD ;
				break ;
				default : 
					cout << "Unknown case in Domain_bispheric_chi_first::base_mult_sin_phi" << endl ;
					abort() ;
				}
		}
	while (index.inc()) ;
	
	*res.base.bases_1d[0] = *so.base.bases_1d[0] ;

	res.base.def = true ;

	res.cf = new Array<double> (so.base.ope_1d(mult_cos_1d, 2, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}

Val_domain Domain_bispheric_chi_first::mult_sin_phi (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base.allocate (nbr_coefs) ;	
	*res.base.bases_1d[2] = *so.base.bases_1d[2] ;
	
	Index index (res.base.bases_1d[1]->get_dimensions()) ;

	// Inversion in theta :
	do {
		switch ((*so.base.bases_1d[1])(index)) {
			case CHEB_EVEN :
				res.base.bases_1d[1]->set(index) = CHEB_ODD ;
				break ;
			case CHEB_ODD :
				res.base.bases_1d[1]->set(index) = CHEB_EVEN ;
				break ;
			case LEG_ODD :
				res.base.bases_1d[1]->set(index) = LEG_EVEN ;
				break ;
			case LEG_EVEN :
				res.base.bases_1d[1]->set(index) = LEG_ODD ;
				break ;
				default : 
					cout << "Unknown case in Domain_bispheric_chi_first::base_mult_cos_phi" << endl ;
					abort() ;
				}
		}
	while (index.inc()) ;
	*res.base.bases_1d[0] = *so.base.bases_1d[0] ;

	res.base.def = true ;

	res.cf = new Array<double> (so.base.ope_1d(mult_sin_1d, 2, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}

Val_domain Domain_bispheric_chi_first::div_chi (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;
	
	res.cf = new Array<double> (so.base.ope_1d(div_x_1d, 1, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}

Val_domain Domain_bispheric_chi_first::div_sin_chi(const Val_domain& so) const {

	// First division in configuration space :
	Val_domain res (so/sin(*p_chi)) ;	

	// Points where division is singular it is done in the coefficient space
	Val_domain copie (so.div_chi()) ;
	copie.coef_i() ;

	Index index (nbr_points) ;
	do {
		if (index(1)==0)
			res.set(index) = copie(index)/chi_max ;
	}
	while (index.inc()) ;


	res.base = copie.base ;
	return res ;
}

Val_domain Domain_bispheric_chi_first::mult_r (const Val_domain& so) const {
    Val_domain res (so*get_radius()) ;
    res.set_base() = so.get_base() ;
    return res ;
}
}
