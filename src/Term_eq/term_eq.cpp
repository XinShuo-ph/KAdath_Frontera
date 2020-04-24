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
#include "tensor_impl.hpp"
#include "space.hpp"
namespace Kadath {
void affecte_one_dom (int, Tensor*, const Tensor*) ;


Term_eq::Term_eq (int dd, const Tensor& vx) : dom{dd}, val_d{nullptr}, der_d{nullptr}, der_t{nullptr}, type_data {TERM_T} {
	val_t = new Tensor{vx, false} ;
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

Term_eq::Term_eq (int dd, const Tensor& vx, const Tensor& dx) :
    dom{dd}, val_d{nullptr}, der_d{nullptr}, val_t{new Tensor{vx, false}}, der_t{new Tensor{dx, false}},
    type_data {TERM_T}
{
	for (int i=0 ; i<val_t->get_n_comp() ; i++) {
		Array<int> id (val_t->indices(i)) ;
		val_t->set(id).set_domain(dom) = vx(id)(dom) ;
		}

	if (vx.is_name_affected()) {
		val_t->set_name_affected() ;
		for (int i=0 ; i<val_t->get_valence() ; i++)
			val_t->set_name_ind(i, vx.get_name_ind()[i]) ;
	}

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

Term_eq::Term_eq (const Term_eq& so) : dom{so.dom}, val_d{nullptr}, der_d{nullptr}, val_t{nullptr}, der_t{nullptr},
						type_data{so.type_data} {

	if (so.val_d!=nullptr)
		val_d = new double{*so.val_d} ;
	if (so.der_d!=nullptr)
		der_d = new double{*so.der_d} ;
	if (so.val_t!=nullptr) {
		val_t = new Tensor {*so.val_t, false} ;
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

	if (so.der_t!=nullptr) {
		der_t = new Tensor {*so.der_t, false} ;
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

Term_eq::Term_eq(Kadath::Term_eq &&so) :
    dom{so.dom},
    val_d{so.val_d},
    der_d{so.der_d},
    val_t{so.val_t},
    der_t{so.der_t},
    type_data{so.type_data}
{
    so.val_d = nullptr;
    so.der_d = nullptr;
    so.val_t = nullptr;
    so.der_t = nullptr;
}
Term_eq& Term_eq::operator=(Term_eq && so)
{
    assert(dom == so.dom && type_data == so.type_data);
    std::swap(val_d,so.val_d);
    std::swap(der_d,so.der_d);
    std::swap(val_t,so.val_t);
    std::swap(der_t,so.der_t);
    return *this;
}

Term_eq::~Term_eq() {

	if (val_d!=nullptr)
		delete val_d ;
	if (der_d!=nullptr)
		delete der_d ;
	if (val_t!=nullptr)
		delete val_t ;
	if (der_t!=nullptr)
		delete der_t ;
}

double Term_eq::get_val_d() const {

	if (type_data!=TERM_D) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (val_d ==nullptr) {
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
	if (der_d ==nullptr) {
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
	if (val_t ==nullptr) {
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
	if (der_t ==nullptr) {
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
		if (so.val_d!=nullptr) {
			if (val_d!=nullptr)
				delete val_d ;
			val_d = new double(*so.val_d) ;
			}
		if (so.der_d!=nullptr) {
			if (der_d!=nullptr)
				delete der_d ;
			der_d = new double(*so.der_d) ;
			}
	}

	if (type_data==TERM_T) {
		
	if (so.val_t!=nullptr) {
	  if (val_t==nullptr)
	    val_t = new Tensor(*so.val_t) ;
	  else
	    affecte_one_dom (dom, val_t, so.val_t) ;
	}
	if (so.der_t!=nullptr) {
	  if (der_t==nullptr)
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
	if (val_d!=nullptr)
		delete val_d ;
	val_d = new double(so) ;
}

void Term_eq::set_der_d (double so) {
	if (type_data!=TERM_D) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (der_d!=nullptr)
		delete der_d ;
	der_d = new double(so) ;
}

void Term_eq::set_val_t (Tensor so) {
	if (type_data!=TERM_T) {
		cerr << "Wrong type of data in Term_eq" << endl ;
		abort() ;
	}
	if (val_t!=nullptr)
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
	if (der_t!=nullptr)
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
			if (der_d!=nullptr)
				delete der_d ;
			der_d = new double(0.) ;
			break ;
		case (TERM_T) :
			assert (val_t!=nullptr) ;
			if (der_t==nullptr)
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
			if (so.val_d !=nullptr)
				flux << "val = " << *so.val_d << endl ;
			else
				flux << "val undefined" << endl ;
			if (so.der_d !=nullptr)
				flux << "der = " << *so.der_d << endl ;
			else
				flux << "der undefined" << endl ;
			break ;
		case (TERM_T) :
			flux << "tensorial data" << endl ;
			if (so.val_t !=nullptr)
				flux << "val = " << *so.val_t << endl ;
			else
				flux << "val undefined" << endl ;
			if (so.der_t !=nullptr)
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
