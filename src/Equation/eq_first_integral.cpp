/*
    Copyright 2018 Philippe Grandclement

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
#include "name_tools.hpp"
namespace Kadath {
Eq_first_integral::Eq_first_integral(const System_of_eqs* syst, const Domain* innerdom, int dommin, int dommax, const char* integ_part, const char* const_part) : Equation(innerdom, dommin, dommax-dommin+2), dom_min(dommin), dom_max(dommax) {


	// First the constant part
	char auxicst[LMAX] ;
	trim_spaces(auxicst, const_part) ;
	parts[0] = syst->give_ope (dom_min, auxicst) ;

	// First all the integral parts
	char auxiinteg[LMAX] ;
	trim_spaces(auxiinteg, integ_part) ;

	for (int d=dom_min ; d<=dom_max; d++)
		parts[d-dom_min+1] = syst->give_ope (d, auxiinteg) ;
}

Eq_first_integral::~Eq_first_integral() {
}


void Eq_first_integral::export_val (int& conte, Term_eq** residus, Array<double>& sec, int& pos_res) const {

	// Standard parts
	assert (residus[conte]->get_type_data()==TERM_T) ;
	// Recover pointer on the space
	const Space& space = residus[conte]->get_val_t().get_space() ;

	// Get the constant part
	Index pori (space.get_domain(dom_min)->get_nbr_points()) ;
	double val_cst_part = residus[conte]->get_val_t()()(dom_min)(pori) ;
	conte ++ ;

	// Recover origin value of first integral
	double val_ori = residus[conte]->get_val_t()()(dom_min)(pori) ;

	// Inner domain
	// Remove value at the origin
	Scalar auximin (space) ;
	auximin.set_domain(dom_min) = residus[conte]->get_val_t()()(dom_min) - val_ori ;
	auximin.set_domain(dom_min).coef_i() ;

	// At the origin put the cst_part (no continuity but does it matter ?)
	Index pos (space.get_domain(dom_min)->get_nbr_points()) ;
	do {
		if (pos(0)==0)
			auximin.set_domain(dom_min).set(pos) = val_cst_part ;
	}
	while (pos.inc()) ;

	space.get_domain(dom_min)->export_tau (auximin, dom_min, 0, sec, pos_res, *n_cond) ;
	conte ++ ;

	// Export all the other domains
	for (int d=dom_min+1 ; d<=dom_max ; d++) {
		Scalar auxi (space) ;
		auxi.set_domain(d) = residus[conte]->get_val_t()()(d) - val_ori ;
		space.get_domain(d)->export_tau (auxi, d, 0, sec, pos_res, *n_cond) ;
		conte ++ ;
	}
}


void Eq_first_integral::export_der (int& conte, Term_eq** residus, Array<double>& sec, int& pos_res) const {
	// Standard parts
	assert (residus[conte]->get_type_data()==TERM_T) ;
	// Recover pointer on the space
	const Space& space = residus[conte]->get_val_t().get_space() ;

	int start = pos_res ;

	// Get the constant part
	Index pori (space.get_domain(dom_min)->get_nbr_points()) ;
	double val_cst_part = residus[conte]->get_der_t()()(dom_min)(pori) ;
	conte ++ ;

	// Recover origin value of first integral
	double val_ori = 0.;
	if(!residus[conte]->get_der_t()()(dom_min).check_if_zero())
		val_ori = residus[conte]->get_der_t()()(dom_min)(pori) ;

	// Inner domain
	// Remove value at the origin
	Scalar auximin (space) ;
	auximin.set_domain(dom_min) = residus[conte]->get_der_t()()(dom_min) - val_ori ;
	auximin.set_domain(dom_min).coef_i() ;

	// At the origin put the cst_part (no continuity but does it matter ?)
	Index pos (space.get_domain(dom_min)->get_nbr_points()) ;
	do {
		if (pos(0)==0)
			auximin.set_domain(dom_min).set(pos) = val_cst_part ;
	}
	while (pos.inc()) ;

	space.get_domain(dom_min)->export_tau (auximin, dom_min, 0, sec, pos_res, *n_cond) ;
	conte ++ ;

	// Export all the other domains
	for (int d=dom_min+1 ; d<=dom_max ; d++) {
		Scalar auxi (space) ;
		auxi.set_domain(d) = residus[conte]->get_der_t()()(d) - val_ori ;
		space.get_domain(d)->export_tau (auxi, d, 0, sec, pos_res, *n_cond) ;
		conte ++ ;
	}
}


Array<int> Eq_first_integral::do_nbr_conditions (const Tensor& tt)  const {

	Array<int> res (dom_max - dom_min+1) ;
	for (int d=dom_min ; d<=dom_max ; d++)
		res.set(d-dom_min) = tt.get_space().get_domain(d)->nbr_conditions (tt, ndom, 0)(0) ;
	return res ;
}


bool Eq_first_integral::take_into_account (int target) const {
	if ((target >= dom_min) && (target<=dom_max))
		return true ;
	else
		return false ;
}

}
