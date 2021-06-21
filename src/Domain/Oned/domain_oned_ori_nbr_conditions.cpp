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
int Domain_oned_ori::nbr_conditions_val_domain (const Val_domain& so, int order) const {
	
	int res = 0 ;
	int basex = (*so.get_base().bases_1d[0])(0) ;
	int max= 0 ;
	switch (basex) {
		case CHEB_EVEN:
			max = nbr_coefs(0) ;
			break ;
		case LEG_EVEN:
			max = nbr_coefs(0) ;
			break ;
		case CHEB_ODD:
			max = nbr_coefs(0)-1 ;
			break ;
		case LEG_ODD:
			max = nbr_coefs(0)-1 ;
			break ;
		default:
			cerr << "Unknow basis in Domain_oned_ori::nbr_conditions_val_domain" << endl ;
			abort() ;
	}

	switch (order) {
		case 0 :
			res = max ;
			break ;
		case 1 :
			res = max-1 ;
			break ;
		case 2 :
			res = max-1 ;
			break ;
		default:
			cerr << "Unknow order in Domain_oned_ori::nbr_conditions_val_domain" << endl ;
			abort() ;
	}
	return res ;
}

Array<int> Domain_oned_ori::nbr_conditions (const Tensor& tt, int dom, int order, int n_cmp, Array<int>**) const {

	int size = (n_cmp==-1) ? tt.get_n_comp() : n_cmp ;
	Array<int> res (size) ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			res.set(0) = nbr_conditions_val_domain (tt()(dom), order) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_oned_ori::nbr_conditions" << endl ;
			break ;
	}
	return res ;
}}
