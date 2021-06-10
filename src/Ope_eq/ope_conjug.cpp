/*
    Copyright 2021 Philippe Grandclement

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

#include "ope_eq.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
namespace Kadath {
Ope_conjug::Ope_conjug (const System_of_eqs* zesys, Ope_eq* target) : Ope_eq(zesys, target->get_dom(), 1) {
	parts[0] = target ;
	
}

Ope_conjug::~Ope_conjug() {
}

Term_eq Ope_conjug::action() const {

	Term_eq target (parts[0]->action()) ;
	
	// Check it is a tensor
	if (target.type_data != TERM_T) {
		cerr << "Ope_conjug only defined with respect for a tensor" << endl ;
		abort() ;
	}

	if (target.val_t->get_n_comp() != 1) {
		cerr << "Ope_conjug only defined with respect to a scalar" << endl ;
		abort() ;
	}

	int m_res = inv_m_quant (target.val_t->get_parameters()) ;
        if (m_res!=0) {
		target.val_t->set_parameters().set_m_quant() = m_res ;
		if (target.der_t!=0x0)
			target.der_t->set_parameters().set_m_quant() = m_res ;
        }
       
	return target ;
}}
