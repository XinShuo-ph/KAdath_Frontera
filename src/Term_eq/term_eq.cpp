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

#include "term_eq.hpp"
#include "scalar.hpp"
#include "space.hpp"
namespace Kadath {
void affecte_one_dom (int, Tensor*, const Tensor*) ;

Term_eq::Term_eq (int dd, int tipe) : dom(dd), val_d(0x0), der_d(0x0), val_t(0x0), der_t(0x0), type_data (tipe) {
	assert ((tipe==TERM_D) || (tipe==TERM_T)) ;
}

Term_eq::Term_eq (int dd, double vx) : dom(dd), der_d(0x0), val_t(0x0), der_t(0x0), type_data (TERM_D) {

	val_d = new double(vx) ;
}

Term_eq::Term_eq (int dd, double vx, double dx) : dom(dd), val_t(0x0), der_t(0x0), type_data (TERM_D) {
	val_d = new double(vx) ;
	der_d = new double(dx) ;
}

Term_eq::Term_eq (int dd, const Tensor& vx) : dom(dd), val_d(0x0), der_d(0x0), der_t(0x0), type_data (TERM_T) {
	val_t = new Tensor(vx, false) ;
	for (int i=0 ; i<val_t->get_n_comp() ; i++) {
		Array<int> id (val_t->indices(i)) ;
		val_t->set(id).set_domain(dom) = vx(id)(dom) ;
		}
	if (vx.is_name_affected()) {
		val_t->set_name_affected() ;
		for (int i=0 ; i<val_t->get_valence() ; i++)
			val_t->set_name_ind(i, vx.get_name_ind()[i]) ;
	}
}

Term_eq::Term_eq (int dd, const Tensor& vx, const Tensor& dx) : dom(dd), val_d(0x0), der_d(0x0), type_data (TERM_T) {
	val_t = new Tensor(vx, false) ;
	for (int i=0 ; i<val_t->get_n_comp() ; i++) {
		Array<int> id (val_t->indices(i)) ;
		val_t->set(id).set_domain(dom) = vx(id)(dom) ;
		}

	if (vx.is_name_affected()) {
		val_t->set_name_affected() ;
		for (int i=0 ; i<val_t->get_valence() ; i++)
			val_t->set_name_ind(i, vx.get_name_ind()[i]) ;
	}

	der_t = new Tensor(dx, false) ;
	for (int i=0 ; i<der_t->get_n_comp() ; i++) {
		Array<int> id (der_t->indices(i)) ;
		der_t->set(id).set_domain(dom) = dx(id)(dom) ;
		}

	if (dx.is_name_affected()) {
		der_t->set_name_affected() ;
		for (int i=0 ; i<der_t->get_valence() ; i++)
			der_t->set_name_ind(i, dx.get_name_ind()[i]) ;
	}
}

Term_eq::Term_eq (const Term_eq& so) : dom(so.dom), val_d(0x0), der_d(0x0), val_t(0x0), der_t(0x0), 
						type_data(so.type_data) {

	if (so.val_d!=0x0)
		val_d = new double(*so.val_d) ;
	if (so.der_d!=0x0)
		der_d = new double(*so.der_d) ;
	if (so.val_t!=0x0) {
		val_t = new Tensor (*so.val_t, false) ;
		for (int i=0 ; i<val_t->get_n_comp() ; i++) {
		Array<int> id (val_t->indices(i)) ;
		val_t->set(id).set_domain(dom) = (*so.val_t)(id)(dom) ;
		}	

	if (so.val_t->is_name_affected()) {
		val_t->set_name_affected() ;
		for (int i=0 ; i<val_t->get_valence() ; i++)
			val_t->set_name_ind(i, so.val_t->get_name_ind()[i]) ;
	}
	}

	if (so.der_t!=0x0) {
		der_t = new Tensor (*so.der_t, false) ;
		for (int i=0 ; i<der_t->get_n_comp() ; i++) {
		Array<int> id (der_t->indices(i)) ;
		der_t->set(id).set_domain(dom) = (*so.der_t)(id)(dom) ;
		}

	if (so.der_t->is_name_affected()) {
		der_t->set_name_affected() ;
		for (int i=0 ; i<der_t->get_valence() ; i++)
			der_t->set_name_ind(i, so.der_t->get_name_ind()[i]) ;
	}
	}
}

Term_eq::~Term_eq() {

	if (val_d!=0x0)
		delete val_d ;
	if (der_d!=0x0)
		delete der_d ;
	if (val_t!=0x0)
		delete val_t ;
	if (der_t!=0x0)
		delete der_t ;
}

double Term_eq::get_val_d() const {

	if (type_data!=TERM_D) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (val_d ==0x0) {
		cerr << "val_d uninitialised in Term_eq" << endl ;
		abort() ;
	}
	return *val_d ;
}

double Term_eq::get_der_d() const {

	if (type_data!=TERM_D) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (der_d ==0x0) {
		cerr << "der_d uninitialised in Term_eq" << endl ;
		abort() ;
	}
	return *der_d ;
}

Tensor Term_eq::get_val_t() const {

	if (type_data!=TERM_T) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (val_t ==0x0) {
		cerr << "val_t uninitialised in Term_eq" << endl ;
		abort() ;
	}
	return *val_t ;
}

Tensor Term_eq::get_der_t() const {

	if (type_data!=TERM_T) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (der_t ==0x0) {
		cerr << "der_t uninitialised in Term_eq" << endl ;
		abort() ;
	}
	return *der_t ;
}

void Term_eq::operator= (const Term_eq& so) {

	assert (dom==so.dom)  ;

	if (type_data!=so.type_data) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}

	if (type_data==TERM_D) {
		if (so.val_d!=0x0) {
			if (val_d!=0x0)
				delete val_d ;
			val_d = new double(*so.val_d) ;
			}
		if (so.der_d!=0x0) {
			if (der_d!=0x0)
				delete der_d ;
			der_d = new double(*so.der_d) ;
			}
	}

	if (type_data==TERM_T) {
		
	if (so.val_t!=0x0) {
	  if (val_t==0x0)
	    val_t = new Tensor(*so.val_t) ;
	  else
	    affecte_one_dom (dom, val_t, so.val_t) ;
	}
	if (so.der_t!=0x0) {
	  if (der_t==0x0)
	    der_t = new Tensor(*so.der_t) ;
	  else
	  affecte_one_dom(dom, der_t, so.der_t) ;
	}
	}
	
}

void Term_eq::set_val_d (double so) {
	if (type_data!=TERM_D) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (val_d!=0x0)
		delete val_d ;
	val_d = new double(so) ;
}

void Term_eq::set_der_d (double so) {
	if (type_data!=TERM_D) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (der_d!=0x0)
		delete der_d ;
	der_d = new double(so) ;
}

void Term_eq::set_val_t (Tensor so) {
	if (type_data!=TERM_T) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (val_t!=0x0)
		delete val_t ;
	val_t = new Tensor(so, false) ;
	for (int i=0 ; i<val_t->get_n_comp() ; i++) {
		Array<int> id (val_t->indices(i)) ;
		val_t->set(id).set_domain(dom) = so(id)(dom) ;
		}
}

void Term_eq::set_der_t (Tensor so) {
	if (type_data!=TERM_T) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (der_t!=0x0)
		delete der_t ;
	der_t = new Tensor(so, false) ;
	for (int i=0 ; i<val_t->get_n_comp() ; i++) {
		Array<int> id (der_t->indices(i)) ;
		der_t->set(id).set_domain(dom) = so(id)(dom) ;
	}
}

void Term_eq::set_der_zero() {

	switch (type_data) {
		case (TERM_D) :
			if (der_d!=0x0)
				delete der_d ;
			der_d = new double(0.) ;
			break ;
		case (TERM_T) :
			assert (val_t!=0x0) ;
			if (der_t==0x0)
				der_t = new Tensor(*val_t, false) ;
			for (int i=0 ; i<der_t->get_n_comp() ; i++)
				der_t->set(der_t->indices(i)).set_domain(dom).set_zero() ;
			break;
		default : 
			cerr << "Wrong type of data in Term_eq" << endl ;
			abort() ;
	}
}

ostream& operator<< (ostream& flux, const Term_eq& so) {
	flux << "Data defined in domain = " << so.dom << endl ;
	switch (so.type_data) {
		case (TERM_D) :
			flux << "double data" << endl ;
			if (so.val_d !=0x0)
				flux << "val = " << *so.val_d << endl ;
			else
				flux << "val undefined" << endl ;
			if (so.der_d !=0x0)
				flux << "der = " << *so.der_d << endl ;
			else
				flux << "der undefined" << endl ;
			break ;
		case (TERM_T) :
			flux << "tensorial data" << endl ;
			if (so.val_t !=0x0)
				flux << "val = " << *so.val_t << endl ;
			else
				flux << "val undefined" << endl ;
			if (so.der_t !=0x0)
				flux << "der = " << *so.der_t << endl ;
			else
				flux << "der undefined" << endl ;
			break ;
		default: 
			cerr << "Unknown data type in Term_eq" << endl ;
			abort() ;
		}
	return flux ;
}}
