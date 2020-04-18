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
#include "tensor_impl.hpp"
namespace Kadath {
Eq_matching_exception::Eq_matching_exception(const Domain* zedom, int dd, int bb, int other_dd, int other_bb, 
								Ope_eq* lhs, Ope_eq* rhs, const Param& par, Ope_eq* exception, int nused, Array<int>** pused) : 
		Equation(zedom, dd, 3, nused, pused), bound(bb),  other_dom(other_dd), other_bound(other_bb), parameters(par)  {
	parts[0] = lhs ;
	parts[1] = rhs ;
	parts[2] = exception ;
}

Eq_matching_exception::~Eq_matching_exception() {
}

void Eq_matching_exception::export_val(int& conte, Term_eq** residus, Array<double>& sec, int& pos_res) const {
	
	assert (residus[conte]->get_type_data()==TERM_T) ;
	assert (residus[conte+1]->get_type_data()==TERM_T) ;
	assert (residus[conte+2]->get_type_data()==TERM_T) ;
	
	int start = pos_res ;
	dom->export_tau_boundary_exception (residus[conte]->get_val_t(), ndom, bound, sec, pos_res, *n_cond, 
						  parameters, 1, residus[conte+2]->get_val_t(), n_cmp_used, p_cmp_used) ;
	Array<double> auxi (pos_res - start) ;
	auxi = 0. ;
	int zero = 0 ;
	residus[conte+1]->get_val_t().get_space().get_domain(other_dom)->export_tau_boundary_exception 
		(residus[conte+1]->get_val_t(), other_dom, other_bound, auxi, zero, *n_cond, parameters, 2, residus[conte+2]->get_val_t(), n_cmp_used, p_cmp_used) ;
	for (int i =start ; i<pos_res ; i++)
		sec.set(i) -= auxi(i-start) ;
	conte +=3 ;
}

void Eq_matching_exception::export_der(int& conte, Term_eq** residus, Array<double>& sec, int& pos_res) const {
	
	assert (residus[conte]->get_type_data()==TERM_T) ;
	assert (residus[conte+1]->get_type_data()==TERM_T) ;
	assert (residus[conte+2]->get_type_data()==TERM_T) ;
	
	int start = pos_res ;
	dom->export_tau_boundary_exception (residus[conte]->get_der_t(), ndom, bound, sec, pos_res, *n_cond, parameters, 3, residus[conte+2]->get_der_t(), 
							n_cmp_used, p_cmp_used) ;
	Array<double> auxi (pos_res - start) ;
	auxi = 0. ;
	int zero = 0 ;
	residus[conte+1]->get_val_t().get_space().get_domain(other_dom)->export_tau_boundary_exception 
		(residus[conte+1]->get_der_t(), other_dom, other_bound, auxi, zero, *n_cond, parameters, 4, residus[conte+2]->get_der_t(), 
															  n_cmp_used, p_cmp_used) ;
	for (int i =start ; i<pos_res ; i++)
		sec.set(i) -= auxi(i-start) ;
	conte +=3 ;
}

Array<int> Eq_matching_exception::do_nbr_conditions (const Tensor& tt)  const {
	return dom->nbr_conditions_boundary (tt, ndom, bound, n_cmp_used, p_cmp_used) ;
}


bool Eq_matching_exception::take_into_account (int target) const {
	if ((target==ndom) || (target==other_dom))
		return true ;
	else
		return false ;
}

}
