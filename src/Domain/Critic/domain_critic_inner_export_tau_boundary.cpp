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
#include "critic.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
void Domain_critic_inner::export_tau_val_domain_boundary (const Val_domain& so, int bound, Array<double>& sec, int& pos_sec, int ncond) const {
	
	if (so.check_if_zero())
		pos_sec += ncond ;
	else {
		so.coef() ;
		int base_t = (*so.get_base().bases_1d[1])(0) ;
		assert ((base_t==COSSIN_EVEN) || (base_t==COSSIN_ODD)) ;

		Index pos_cf (nbr_coefs) ;
		// Loop on t
		int max = (base_t==COSSIN_EVEN) ? nbr_coefs(1)-1 : nbr_coefs(1)-2 ;
		for (int j=0 ; j<max ; j++)
			if ((j!=1) || (base_t!=COSSIN_EVEN)) {
				pos_cf.set(1) = j ;
				sec.set(pos_sec) = val_boundary (bound, so, pos_cf) ;
				pos_sec ++ ;
			}
	}
}

void Domain_critic_inner::export_tau_boundary (const Tensor& tt, int dom, int bound, Array<double>& res, int& pos_res, const Array<int>& ncond,
										int, Array<int>**) const {

	// Check boundary
	if ((bound!=INNER_BC) && (bound!=OUTER_BC)) {
		cerr << "Unknown boundary in Domain_critic_inner::export_tau_boundary" << endl ;
		abort() ;
	}

	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			export_tau_val_domain_boundary (tt()(dom), bound, res, pos_res, ncond(0)) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_critic_inner::export_tau_boundary" << endl ;
			break ;
	}
}}
