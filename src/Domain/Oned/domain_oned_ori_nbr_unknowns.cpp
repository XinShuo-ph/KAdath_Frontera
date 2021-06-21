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

#include "headcpp.hpp"
#include "oned.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
int Domain_oned_ori::nbr_unknowns_val_domain (const Val_domain& so) const {
	
	int res = 0 ;
	int basex = (*so.get_base().bases_1d[0])(0) ;
	switch (basex) {
		case CHEB_EVEN:
			res = nbr_coefs(0) ;
			break ;
		case LEG_EVEN:
			res = nbr_coefs(0) ;
			break ;
		case CHEB_ODD:
			res = nbr_coefs(0)-1 ;
			break ;
		case LEG_ODD:
			res = nbr_coefs(0)-1 ;
			break ;
		default:
			cerr << "Unknowbasis in Domain_oned_ori::nbr_unknowns_val_domain" << endl ;
			abort() ;
	}
	return res ;
}

int Domain_oned_ori::nbr_unknowns (const Tensor& tt, int dom) const {

	// Check right domain
	assert (tt.get_space().get_domain(dom)==this) ;

	int res = 0 ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			res += nbr_unknowns_val_domain (tt()(dom)) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_oned_ori::nbr_unknowns" << endl ;
			break ;
	}
	return res ;
}}
