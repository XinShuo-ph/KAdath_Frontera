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
#include "polar.hpp"
#include "point.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
int Domain_polar_nucleus::nbr_conditions_val_domain_boundary (const Val_domain& so, int mquant) const {
	
	int res = 0 ;
	
	for (int j=0 ; j<nbr_coefs(1) ; j++) {
		bool indic = true ;
		// Get base in theta :
		int baset = (*so.get_base().bases_1d[1])(0) ;
		switch (baset) {
					case COS_EVEN:
						if ((j==0) && (mquant!=0))
							indic = false ;
						break ;
					case COS_ODD:
						if ((j==nbr_coefs(1)-1) || ((j==0) && (mquant!=0)))
							indic = false ;
						break ;
					case SIN_EVEN:
						if (((j==1) && (mquant>1)) ||(j==0) || (j==nbr_coefs(1)-1)) 
							indic = false  ;
						break ;
					case SIN_ODD:
						if (((j==0) && (mquant>1)) ||  (j==nbr_coefs(1)-1))
							indic = false ;
						break ;
					default:
						cerr << "Unknow theta basis in Domain_polar_nucleus::nbr_conditions_val_boundary" << endl ;
						abort() ;
		}

		if (indic)
			res ++ ;
	}
	return res ;
}

Array<int> Domain_polar_nucleus::nbr_conditions_boundary (const Tensor& tt, int dom, int bound, int n_cmp, Array<int>**) const {

	// Check boundary
	if (bound!=OUTER_BC) {
		cerr << "Unknown boundary in Domain_polar_nucleus::nbr_conditions_boundary" << endl ;
		abort() ;
	}

	int size = (n_cmp==-1) ? tt.get_n_comp() : n_cmp ;
	Array<int> res (size) ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			if (!tt.is_m_quant_affected())
			  res.set(0) = nbr_conditions_val_domain_boundary (tt()(dom), 0) ;
			else
			  res.set(0) = nbr_conditions_val_domain_boundary (tt()(dom), tt.get_parameters()->get_m_quant()) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_polar_nucleus::nbr_conditions_boundary" << endl ;
			break ;
	}
	return res ;
}}
