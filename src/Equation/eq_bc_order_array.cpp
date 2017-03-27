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

#include "system_of_eqs.hpp"
#include "ope_eq.hpp"
#include "scalar.hpp"
namespace Kadath {
Eq_bc_order_array::Eq_bc_order_array(const Domain* zedom, int dd, int bb, const Array<int>& ord, Ope_eq* so, int nused, Array<int>** pused) : Equation(zedom, dd, 1, nused, pused), bound(bb), order(ord) {
	parts[0] = so ;
}

Eq_bc_order_array::~Eq_bc_order_array() {
}

void Eq_bc_order_array::export_val(int& conte, Term_eq** residus, Array<double>& sec, int& pos_res) const {
	
	assert (residus[conte]->get_type_data()==TERM_T) ;
	dom->export_tau_boundary_array (residus[conte]->get_val_t(), ndom, bound, order, sec, pos_res, *n_cond, n_cmp_used, p_cmp_used) ;
	conte ++ ;
}

void Eq_bc_order_array::export_der(int& conte, Term_eq** residus, Array<double>& sec, int& pos_res) const {
	
	assert (residus[conte]->get_type_data()==TERM_T) ;
	dom->export_tau_boundary_array (residus[conte]->get_der_t(), ndom, bound, order, sec, pos_res, *n_cond, n_cmp_used, p_cmp_used) ;
	conte ++ ;
}

Array<int> Eq_bc_order_array::do_nbr_conditions (const Tensor& tt)  const {
	return dom->nbr_conditions_boundary_array (tt, ndom, bound, order, n_cmp_used, p_cmp_used) ;
}

bool Eq_bc_order_array::take_into_account (int target) const {
	if (target==ndom)
		return true ;
	else
		return false ;
}
}
