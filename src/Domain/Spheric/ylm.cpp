/*
    Copyright 2022 Philippe Grandclement

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

#include "headcpp.hpp"

#include "utilities.hpp"
#include "spheric.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "tensor.hpp"
#include "term_eq.hpp"
namespace Kadath {


void Domain_shell::ylm_leg_even (Array<double>& cf) const {
	
	static int np_done = 0 ;
	static int nt_done = 0 ;
	static Array<double>* p_mat = 0x0 ;
	
	// Check whether the passage matrix must be computed or not
	int np = nbr_coefs(2)  ;
	int nt = nbr_coefs(1) ;
	int nr = nbr_coefs(0) ;
	
	if ((np!=np_done) || (nt!=nt_done)) {
		np_done = np ;
		nt_done = nt ;
		
		if (p_mat != 0x0)
			delete p_mat ;
		p_mat = new Array<double> (mat_leg_even(nt,np-2)) ;
		
	}
	Array<double> auxi (cf) ;
	cf = 0 ;

	for (int k=0 ; k<np ; k++) // Loop on phi
    	 	if (k!=1) {
    			int m = (k%2==0) ? k/2 : (k-1)/2 ;
    				for (int l=0 ; l<nt ; l++) // Loop on theta
    					for (int j=0 ; j<nt ; j++) // Matrix summation loop 
    						for (int i=0 ; i<nr ; i++) // for all the r
 					   		cf.set(i,l,k) += (*p_mat)(m,l,j)*auxi(i,j,k) ;  
    }
}

void Domain_shell::ylm_inv_leg_even (Array<double>& cf) const {
	
	static int np_done = 0 ;
	static int nt_done = 0 ;
	static Array<double>* p_mat = 0x0 ;
	
	// Check whether the passage matrix must be computed or not
	int np = nbr_coefs(2) ;
	int nt = nbr_coefs(1) ;
	int nr = nbr_coefs(0) ;
	
	if ((np!=np_done) || (nt!=nt_done)) {
		np_done = np ;
		nt_done = nt ;
		
		if (p_mat != 0x0)
			delete p_mat ;
		p_mat = new Array<double> (mat_inv_leg_even(nt,np-2)) ;
	}
	Array<double> auxi (cf) ;
	cf = 0 ;
	
	for (int k=0 ; k<np ; k++) // Loop on phi
    	 	if (k!=1) {
    			int m = (k%2==0) ? k/2 : (k-1)/2 ;
    				for (int l=0 ; l<nt ; l++) // Loop on theta
    					for (int j=0 ; j<nt ; j++) // Matrix summation loop 
    						for (int i=0 ; i<nr ; i++) // for all the r
 					   		cf.set(i,j,k) += (*p_mat)(m,l,j)*auxi(i,l,k) ;  
    }
}

void Domain_shell::ylm_leg_odd (Array<double>& cf) const {
	
	static int np_done = 0 ;
	static int nt_done = 0 ;
	static Array<double>* p_mat = 0x0 ;
	
	// Check whether the passage matrix must be computed or not
	int np = nbr_coefs(2)  ;
	int nt = nbr_coefs(1) ;
	int nr = nbr_coefs(0) ;
	
	if ((np!=np_done) || (nt!=nt_done)) {
		np_done = np ;
		nt_done = nt ;
		
		if (p_mat != 0x0)
			delete p_mat ;
		p_mat = new Array<double> (mat_leg_odd(nt,np-2)) ;
		
	}
	Array<double> auxi (cf) ;
	cf = 0 ;

	for (int k=0 ; k<np ; k++) // Loop on phi
    	 	if (k!=1) {
    			int m = (k%2==0) ? k/2 : (k-1)/2 ;
    				for (int l=0 ; l<nt ; l++) // Loop on theta
    					for (int j=0 ; j<nt ; j++) // Matrix summation loop 
    						for (int i=0 ; i<nr ; i++) // for all the r
 					   		cf.set(i,l,k) += (*p_mat)(m,l,j)*auxi(i,j,k) ;  
    }
}

void Domain_shell::ylm_inv_leg_odd (Array<double>& cf) const {
	
	static int np_done = 0 ;
	static int nt_done = 0 ;
	static Array<double>* p_mat = 0x0 ;
	
	// Check whether the passage matrix must be computed or not
	int np = nbr_coefs(2) ;
	int nt = nbr_coefs(1) ;
	int nr = nbr_coefs(0) ;
	
	if ((np!=np_done) || (nt!=nt_done)) {
		np_done = np ;
		nt_done = nt ;
		
		if (p_mat != 0x0)
			delete p_mat ;
		p_mat = new Array<double> (mat_inv_leg_odd(nt,np-2)) ;
	}
	Array<double> auxi (cf) ;
	cf = 0 ;
	
	for (int k=0 ; k<np ; k++) // Loop on phi
    	 	if (k!=1) {
    			int m = (k%2==0) ? k/2 : (k-1)/2 ;
    				for (int l=0 ; l<nt ; l++) // Loop on theta
    					for (int j=0 ; j<nt ; j++) // Matrix summation loop 
    						for (int i=0 ; i<nr ; i++) // for all the r
 					   		cf.set(i,j,k) += (*p_mat)(m,l,j)*auxi(i,l,k) ;  
    }
}


void Term_eq::ylm () {

	// Check if so is a tensor
	if (type_data != TERM_T) {
		cerr << "Term_eq::ylm only defined for a tensor" << endl ;
		abort() ;
	}

	//Check tensorial basis
	double valence = val_t->get_valence() ;
	if (valence >0)
		if (val_t->get_basis().get_basis(dom) != CARTESIAN_BASIS) {
			cerr << "Term_eq::ylm only defined for a tensor" << endl ;
			abort() ;
	}
	
	bool doder = (der_t==0x0) ? false : true ;

	if (valence==0) {
		val_t->set().at(dom).coef() ;
		// Kill config space
		val_t->set().set_domain(dom).set_in_coef() ;
		// Assumes a real scalar field
		val_t->get_space().get_domain(dom)->ylm_leg_even (*val_t->set().set_domain(dom).cf) ;
		
		if (doder) {
			der_t->set().at(dom).coef() ;
			// Kill config space
			der_t->set().set_domain(dom).set_in_coef() ;
			// Assumes a real scalar field
			der_t->get_space().get_domain(dom)->ylm_leg_even (*val_t->set().set_domain(dom).cf) ;
		}
	}
	
	else {
		Index pos (*val_t) ;
		do {
			// Odd or even
			int sym = 1 ;
			for (int cmp=0 ; cmp<valence ; cmp++)
				if (pos(cmp)==2)
					sym *= -1 ;
					
			val_t->set(pos).at(dom).coef() ;
			val_t->set(pos).set_domain(dom).set_in_coef() ;
			if (sym==1) 
				val_t->get_space().get_domain(dom)->ylm_leg_even (*val_t->set(pos).set_domain(dom).cf) ;
			else 
				val_t->get_space().get_domain(dom)->ylm_leg_odd (*val_t->set(pos).set_domain(dom).cf) ;
				
			if (doder) {
				der_t->set(pos).at(dom).coef() ;
				der_t->set(pos).set_domain(dom).set_in_coef() ;
				if (sym==1) 
					der_t->get_space().get_domain(dom)->ylm_leg_even (*val_t->set(pos).set_domain(dom).cf) ;
				else 
					der_t->get_space().get_domain(dom)->ylm_leg_odd (*val_t->set(pos).set_domain(dom).cf) ;
			}
		}
		while (pos.inc()) ;
	}
}


void Term_eq::ylm_i () {

	// Check if so is a tensor
	if (type_data != TERM_T) {
		cerr << "Term_eq::ylm_i only defined for a tensor" << endl ;
		abort() ;
	}

	//Check tensorial basis
	double valence = val_t->get_valence() ;
	if (valence >0)
		if (val_t->get_basis().get_basis(dom) != CARTESIAN_BASIS) {
			cerr << "Term_eq::ylm_i only defined for a tensor" << endl ;
			abort() ;
	}

	bool doder = (der_t==0x0) ? false : true ;

	if (valence==0) {
		// Assumes a real scalar field
		val_t->get_space().get_domain(dom)->ylm_inv_leg_even (*val_t->set().set_domain(dom).cf) ;
		
		if (doder) 
			der_t->get_space().get_domain(dom)->ylm_inv_leg_even (*val_t->set().set_domain(dom).cf) ;
	}
	
	else {
		Index pos (*val_t) ;
		do {
			// Odd or even
			int sym = 1 ;
			for (int cmp=0 ; cmp<valence ; cmp++)
				if (pos(cmp)==2)
					sym *= -1 ;
			
			if (sym==1) 
				val_t->get_space().get_domain(dom)->ylm_inv_leg_even (*val_t->set(pos).set_domain(dom).cf) ;
			else 
				val_t->get_space().get_domain(dom)->ylm_inv_leg_odd (*val_t->set(pos).set_domain(dom).cf) ;
				
			if (doder)
				if (sym==1) 
					der_t->get_space().get_domain(dom)->ylm_inv_leg_even (*val_t->set(pos).set_domain(dom).cf) ;
				else 
					der_t->get_space().get_domain(dom)->ylm_inv_leg_odd (*val_t->set(pos).set_domain(dom).cf) ;
		}
		while (pos.inc()) ;
	}
}
}
