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
#include "tensor_impl.hpp"
namespace Kadath {
Eq_int::Eq_int (int np) : n_ope(np) {
	parts = new Ope_eq*[n_ope] ;
	for (int i=0 ; i<n_ope ; i++)
		parts[i] = 0x0 ;
}

Eq_int::~Eq_int() {
	for (int i=0 ; i<n_ope ; i++)
		if (parts[i] !=0x0)
			delete parts[i] ;
	delete [] parts ;
}

double Eq_int::get_val() const {
	double res = 0. ;
	for (int i=0 ; i<n_ope ; i++) {
	  res += parts[i]->action().get_val_d() ;
	}
	return res ;
}

double Eq_int::get_der() const {
	double res = 0. ;
	for (int i=0 ; i<n_ope ; i++) 
		res += parts[i]->action().get_der_d() ;
	return res ;
}

void Eq_int::set_part (int pos, Ope_eq* so) {
	assert ((pos>=0) && (pos<n_ope)) ;
	parts[pos] = so ;
}}
