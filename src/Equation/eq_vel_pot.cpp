/*
    Copyright 2019 Philippe Grandclement

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
Eq_vel_pot::Eq_vel_pot(const Domain* zedom, int dd, int ord, Ope_eq* so, Ope_eq* constant) : Equation(zedom, dd, 2), order(ord) {
	parts[0] = so ;
	parts[1] = constant ;
}

Eq_vel_pot::~Eq_vel_pot() {
}

void Eq_vel_pot::export_val (int& conte, Term_eq** residus, Array<double>& sec, int& pos_res) const {
	
	assert (residus[conte]->get_type_data()==TERM_T) ;
	assert (residus[conte]->get_val_t().get_valence()==0) ; //only defined for a scalar field so far
	
	assert (residus[conte+1]->get_type_data()==TERM_T) ;
	assert (residus[conte+1]->get_val_t().get_valence()==0) ; 

	int old_pos = pos_res ;

	dom->export_tau (residus[conte]->get_val_t(), ndom, order, sec, pos_res, *n_cond) ;
	
	// Get the first coef (a bit long probably but anyway...)
	Array<double> auxi (pos_res - old_pos) ;
	auxi = 0. ;
	int zero = 0 ;
	dom->export_tau(residus[conte+1]->get_val_t(), ndom, order, auxi, zero, *n_cond) ;
	
	// Put the coef in the right place 
	sec.set(old_pos) = auxi(0) ;

	conte += 2 ;
}

void Eq_vel_pot::export_der (int& conte, Term_eq** residus, Array<double>& sec, int& pos_res) const {
	
	assert (residus[conte]->get_type_data()==TERM_T) ;
	assert (residus[conte]->get_der_t().get_valence()==0) ; //only defined for a scalar field so far
	
	assert (residus[conte+1]->get_type_data()==TERM_T) ;
	assert (residus[conte+1]->get_der_t().get_valence()==0) ; 

	int old_pos = pos_res ;

	dom->export_tau (residus[conte]->get_der_t(), ndom, order, sec, pos_res, *n_cond) ;
	
	// Get the first coef (a bit long probably but anyway...)
	Array<double> auxi (pos_res - old_pos) ;
	auxi = 0. ;
	int zero = 0 ;
	dom->export_tau(residus[conte+1]->get_der_t(), ndom, order, auxi, zero, *n_cond) ;
	
	// Put the coef in the right place 
	sec.set(old_pos) = auxi(0) ;

	
	conte +=2  ;
}

Array<int> Eq_vel_pot::do_nbr_conditions (const Tensor& tt)  const {
	return dom->nbr_conditions (tt, ndom, order) ;
}

bool Eq_vel_pot::take_into_account (int target) const {
	if (target == ndom)
		return true ;
	else 
		return false ;
}

}
