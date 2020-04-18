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
#include "array_math.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "tensor.hpp"

namespace Kadath {
int Domain_oned_inf::nbr_conditions_val_domain_boundary (const Val_domain&) const {
	// Pretty simple in the 1d case...
	int res = 1 ;
	return res ;
}

Array<int> Domain_oned_inf::nbr_conditions_boundary (const Tensor& tt, int dom, int bound, int n_cmp, Array<int>**) const {

	// Check boundary
	if ((bound!=OUTER_BC) && (bound!=INNER_BC)) {
		cerr << "Unknown boundary in Domain_oned_inf::nbr_conditions_boundary" << endl ;
		abort() ;
	}

	int size = (n_cmp==-1) ? tt.get_n_comp() : n_cmp ;
	Array<int> res (size) ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			res.set(0) = nbr_conditions_val_domain_boundary (tt()(dom)) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_oned_inf::nbr_conditions_boundary" << endl ;
			break ;
	}
	return res ;
}}
