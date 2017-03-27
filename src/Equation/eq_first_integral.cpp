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
#include "term_eq.hpp"
#include "scalar.hpp"
namespace Kadath {
Eq_first_integral::Eq_first_integral(const Domain* zedom, int dd, Ope_eq* so, Ope_eq* ori, int nused, Array<int>** pused) : Equation(zedom, dd, 2, nused, pused) {
	parts[0] = so ;
	parts[1] = ori ;
}

Eq_first_integral::~Eq_first_integral() {
}

void Eq_first_integral::export_val (int& conte, Term_eq** residus, Array<double>& sec, int& pos_res) const {
	
	// Standard part 
	assert (residus[conte]->get_type_data()==TERM_T) ;
	int start = pos_res ;
	dom->export_tau (residus[conte]->get_val_t(), ndom, 0, sec, pos_res, *n_cond, n_cmp_used, p_cmp_used) ;
	conte ++ ;
	
	// Integral part :
	assert (residus[conte]->get_type_data()==TERM_D) ;
	sec.set(start) -= residus[conte]->get_val_d() ;
	conte ++ ;
}

void Eq_first_integral::export_der (int& conte, Term_eq** residus, Array<double>& sec, int& pos_res) const {
	
	assert (residus[conte]->get_type_data()==TERM_T) ;
	int start = pos_res ;
	dom->export_tau (residus[conte]->get_der_t(), ndom, 0, sec, pos_res, *n_cond, n_cmp_used, p_cmp_used) ;
	conte ++ ;
	
	// Integral part :
	assert (residus[conte]->get_type_data()==TERM_D) ;
	sec.set(start) -= residus[conte]->get_der_d() ;
	conte ++ ;
}

Array<int> Eq_first_integral::do_nbr_conditions (const Tensor& tt)  const {
	return dom->nbr_conditions (tt, ndom, 0, n_cmp_used, p_cmp_used) ;
}

bool Eq_first_integral::take_into_account (int target) const {
	if ((target == ndom) || (target==0))
		return true ;
	else 
		return false ;
}

}
