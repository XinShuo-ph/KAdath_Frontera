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

#include "ope_eq.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "headcpp.hpp"

namespace Kadath {
Ope_sin::Ope_sin (const System_of_eqs* zesys, Ope_eq* target) : Ope_eq(zesys, target->get_dom(), 1) {
	parts[0] = target ;
}

Ope_sin::~Ope_sin() {
}

Term_eq Ope_sin::action() const {

	Term_eq target (parts[0]->action()) ;
	switch (target.type_data) {
	    case TERM_T : {
 
	    // Check it is a scalar
	    if (target.val_t->get_valence() != 0) {
		cerr << "Ope_sin only defined with respect for a scalar" << endl ;
		abort() ;
	    }
	    
	    // The value	
	Tensor resval (*target.val_t, false) ;

	for (int i=0 ; i<target.val_t->get_n_comp() ; i++) {
		Array<int> ind (target.val_t->indices(i)) ;
		Val_domain value ((*target.val_t)(ind)(dom)) ;
		resval.set(ind).set_domain(dom) = sin(value) ;
	}

	if (target.der_t!=0x0) {
		Tensor resder (*target.der_t, false) ;
		for (int i=0 ; i<target.der_t->get_n_comp() ; i++) {
			Array<int> ind (target.der_t->indices(i)) ;
			Val_domain valder ((*target.der_t)(ind)(dom)) ;
			Val_domain value ((*target.val_t)(ind)(dom)) ;
			resder.set(ind).set_domain(dom) = valder * cos(value) ;
			}
		Term_eq res (dom, resval, resder) ;
		return res ;
	}
	else {
		Term_eq res (dom, resval) ;
		return res ;
	}
	}
	break ;
	case TERM_D : {
	      if (target.der_d==0x0) {
		  Term_eq res (dom, sin(*target.val_d)) ;
		  return res ;
	      }
	      else {
		  Term_eq res (dom, sin(*target.val_d), *target.der_d * cos(*target.val_d)) ;
		  return res ;
	      }
	  }
	  break ;
	default : {
	  cerr << "Unknown storage in Term_eq..." << endl ;
	  abort() ;
	}
      }
    cerr << "Warning should not be here in Ope_sin::action..." << endl ;
    abort() ;
}
}
