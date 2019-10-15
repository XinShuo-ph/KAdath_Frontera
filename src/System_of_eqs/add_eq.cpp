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
#include "name_tools.hpp"
#include "list_comp.hpp"

namespace Kadath {
 void System_of_eqs::add_eq_inside (int dom, const char* nom, int n_cmp, Array<int>** p_cmp) {
    // Is it written like =0 ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;
		
		// no lhs :
		if (indic)
			eq[neq] = new Eq_inside(espace.get_domain(dom), dom, give_ope(dom, p1), n_cmp, p_cmp) ;
		
		else 
			eq[neq] = new Eq_inside(espace.get_domain(dom), dom, 
						new Ope_sub(this, give_ope(dom, p1), give_ope(dom, p2)), n_cmp, p_cmp) ;
		neq ++ ;
	}
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_inside (int dom, const char* nom, const List_comp& list) {
	add_eq_inside (dom, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_order (int dom, int order, const char* nom, int n_cmp, Array<int>** p_cmp) {
    // Is it written like =0 ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;
		
		// no lhs :
		if (indic)
			eq[neq] = new Eq_order(espace.get_domain(dom), dom, order, give_ope(dom, p1), n_cmp, p_cmp) ;
		
		else 
			eq[neq] = new Eq_order(espace.get_domain(dom), dom, order,
						new Ope_sub(this, give_ope(dom, p1), give_ope(dom, p2)), n_cmp, p_cmp) ;
		neq ++ ;
	}
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_vel_pot (int dom, int order, const char* nom, const char* const_part) {
    // Is it written like =0 ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic1 = is_ope_bin(nom, p1, p2, '=') ;
	
	char p3[LMAX] ;
	char p4[LMAX] ;
	bool indic2 = is_ope_bin(const_part, p3, p4, '=') ;
	

	if ((!indic1) || (!indic2)) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Verif lhs1 = 0 ?
		indic1 = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;

		indic2 = ((p4[0]=='0') && (p4[1]==' ') && (p4[2]=='\0')) ?
			true : false ;
		
		// no lhs :
		if ((indic1) && (indic2))
			eq[neq] = new Eq_vel_pot(espace.get_domain(dom), dom, order, give_ope(dom, p1), give_ope(dom,p3)) ;
		// lhs in 1
		if ((!indic1) && (indic2))
			eq[neq] = new Eq_vel_pot(espace.get_domain(dom), dom, order, new Ope_sub(this, give_ope(dom, p1), give_ope(dom, p2)), give_ope(dom,p3)) ;
		// lhs in 2
		if ((indic1) && (!indic2))
			eq[neq] = new Eq_vel_pot(espace.get_domain(dom), dom, order, give_ope(dom,p1), new Ope_sub(this, give_ope(dom, p3), give_ope(dom, p4))) ;
		// both lhs 
		if ((!indic1) && (!indic2))
			eq[neq] = new Eq_vel_pot(espace.get_domain(dom), dom, order, new Ope_sub(this, give_ope(dom, p1), give_ope(dom, p2)), new Ope_sub(this, give_ope(dom, p3), give_ope(dom, p4))) ;

		neq ++ ;
	}
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_order (int dom, int order, const char* nom, const List_comp& list) {
	add_eq_order (dom, order, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_bc (int dom, int bound, const char* nom, int n_cmp, Array<int>** p_cmp) {
    // Is it written like =0 ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for boundary conditions" << endl ;
		abort() ;
	}
	else {
		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;
			
		// no lhs :
		if (indic)
			eq[neq] = new Eq_bc(espace.get_domain(dom), dom, bound, give_ope(dom, p1, bound), n_cmp, p_cmp) ;
		else
			eq[neq] = new Eq_bc(espace.get_domain(dom), dom, 
				bound, new Ope_sub(this, give_ope(dom, p1, bound), give_ope(dom, p2, bound)), n_cmp, p_cmp) ;

		neq ++ ;
	}	
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_bc (int dom, int bound, const char* nom, const List_comp& list) {
	add_eq_bc (dom, bound, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_matching (int dom, int bound, const char* nom, int n_cmp, Array<int>** p_cmp) {
      int other_dom ;
	int other_bound ;
	espace.get_domain(dom)->find_other_dom (dom, bound, other_dom, other_bound) ;
	assert (other_dom>=dom_min) ;
	assert (other_dom<=dom_max) ;

	// Is it written with =  ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;

	if (!indic) {
		char auxi[LMAX] ;
		trim_spaces(auxi, nom) ;
		// Version without =
		eq[neq] = new Eq_matching(espace.get_domain(dom), dom, bound, other_dom, other_bound, 
			give_ope(dom, auxi, bound), give_ope(other_dom, auxi, other_bound), n_cmp, p_cmp) ;
		neq++ ;
	}
	else {
		// Version with = 
		eq[neq] = new Eq_matching(espace.get_domain(dom), dom, bound, other_dom, other_bound, 
			give_ope(dom, p1, bound), give_ope(other_dom, p2, other_bound), n_cmp, p_cmp) ;
		neq++ ;
	}	
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_matching (int dom, int bound, const char* nom, const List_comp& list) {
	add_eq_matching (dom, bound, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_matching_exception (int dom, int bound, const char* nom, const Param& par, const char* nom_exception, int n_cmp, Array<int>** p_cmp) {
      int other_dom ;
	int other_bound ;
	espace.get_domain(dom)->find_other_dom (dom, bound, other_dom, other_bound) ;
	assert (other_dom>=dom_min) ;
	assert (other_dom<=dom_max) ;

	// Is it written with =  ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;

	if (!indic) {
		char auxi[LMAX] ;
		trim_spaces(auxi, nom) ;
		
		char auxi_exception[LMAX] ;
		trim_spaces(auxi_exception, nom_exception) ;
		// Version without =
		eq[neq] = new Eq_matching_exception(espace.get_domain(dom), dom, bound, other_dom, other_bound, 
			give_ope(dom, auxi, bound), give_ope(other_dom, auxi, other_bound), par, give_ope(dom, auxi_exception, bound), n_cmp, p_cmp) ;
		neq++ ;
	}
	else {
		char auxi_exception[LMAX] ;
		trim_spaces(auxi_exception, nom_exception) ;
	  
		// Version with = 
		eq[neq] = new Eq_matching_exception(espace.get_domain(dom), dom, bound, other_dom, other_bound, 
			give_ope(dom, p1, bound), give_ope(other_dom, p2, other_bound), par, give_ope(dom, auxi_exception, bound) , n_cmp, p_cmp) ;
		neq++ ;
	}	
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_matching_exception (int dom, int bound, const char* nom, const Param& par, const char* nom_exception, const List_comp& list) {
	add_eq_matching_exception (dom, bound, nom, par, nom_exception, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_matching_one_side (int dom, int bound, const char* nom, int n_cmp, Array<int>** p_cmp) {
      int other_dom ;
	int other_bound ;
	espace.get_domain(dom)->find_other_dom (dom, bound, other_dom, other_bound) ;
	assert (other_dom>=dom_min) ;
	assert (other_dom<=dom_max) ;

	// Is it written with =  ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;

	if (!indic) {
		char auxi[LMAX] ;
		trim_spaces(auxi, nom) ;
		// Version without =
		eq[neq] = new Eq_matching_one_side(espace.get_domain(dom), dom, bound, other_dom, other_bound, 
			give_ope(dom, auxi, bound), give_ope(other_dom, auxi, other_bound), n_cmp, p_cmp) ;
		neq++ ;
	}
	else {
		// Version with = 
		eq[neq] = new Eq_matching_one_side(espace.get_domain(dom), dom, bound, other_dom, other_bound, 
			give_ope(dom, p1, bound), give_ope(other_dom, p2, other_bound), n_cmp, p_cmp) ;
		neq++ ;
	}	
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_matching_one_side (int dom, int bound, const char* nom,  const List_comp& list) {
	add_eq_matching_one_side (dom, bound, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_matching_non_std (int dom, int bound, const char* nom, int n_cmp, Array<int>** p_cmp) {

	// First get the number, the indices and associated boundaries of the other domains (member of espace) :
	Array<int> other_props (espace.get_indices_matching_non_std (dom, bound)) ;

	// The equation
	eq[neq] = new Eq_matching_non_std(espace.get_domain(dom), dom, bound, other_props, n_cmp, p_cmp) ;

	// Affectation of the operator in each concerned domain :
	// Current one :bool indic = is_ope_bin(nom, p1, p2, '=') ;
	// Is it written with =  ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;

	if (!indic) {
		char auxi[LMAX] ;
		trim_spaces(auxi, nom) ;
		// Version without =
		eq[neq]->parts[0] = give_ope (dom, auxi, bound) ;
		// The associated ones :
		for (int i=0 ; i<eq[neq]->n_ope-1 ; i++)
			eq[neq]->parts[i+1] = give_ope (other_props(0,i), auxi, other_props(1,i)) ;
		neq++ ;
	}
	else {
		// Version without =
		eq[neq]->parts[0] = give_ope (dom, p1, bound) ;
		// The associated ones :
		for (int i=0 ; i<eq[neq]->n_ope-1 ; i++)
			eq[neq]->parts[i+1] = give_ope (other_props(0,i), p2, other_props(1,i)) ;
		neq++ ;
	}
	
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_matching_non_std (int dom, int bound, const char* nom, const List_comp& list) {
	add_eq_matching_non_std (dom, bound, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_matching_import (int dom, int bound, const char* nom, int n_cmp, Array<int>** p_cmp) {

	// First get the number, the indices and associated boundaries of the other domains (member of espace) :
	Array<int> others (espace.get_indices_matching_non_std (dom, bound)) ;

	// Is it written with =  ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;

	if (!indic) {
		char auxi[LMAX] ;
		trim_spaces(auxi, nom) ;
		//Version without = ; assumes p2 = import(p1)
		eq[neq] = new Eq_matching_import (espace.get_domain(dom), dom, bound, new Ope_sub (this, give_ope(dom, auxi, bound), 
			new Ope_import(this, dom, bound, auxi)) , others, n_cmp, p_cmp) ;
		neq ++ ;
	}
	else {
		// Version with = 
		eq[neq] = new Eq_matching_import (espace.get_domain(dom), dom, bound, new Ope_sub (this, give_ope(dom, p1, bound), give_ope(dom, p2, bound)), others, n_cmp, p_cmp) ;
		neq++ ;
	}	
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_matching_import (int dom, int bound, const char* nom, const List_comp& list) {
	add_eq_matching_import (dom, bound, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_full (int dom, const char* nom, int n_cmp, Array<int>** p_cmp) {

	// Is it written like =0 ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;
		
		// no lhs :
		if (indic)
			eq[neq] = new Eq_full(espace.get_domain(dom), dom, give_ope(dom, p1), n_cmp, p_cmp) ;
		
		else 
			{
			eq[neq] = new Eq_full(espace.get_domain(dom), dom, 
						new Ope_sub(this, give_ope(dom, p1), give_ope(dom, p2)), n_cmp, p_cmp) ;
		}
		neq ++ ;
	}
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_full (int dom, const char* nom, const List_comp& list) {
	add_eq_full (dom, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_one_side (int dom, const char* nom, int n_cmp, Array<int>** p_cmp) {

	// Is it written like =0 ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;
		
		// no lhs :
		if (indic)
			eq[neq] = new Eq_one_side(espace.get_domain(dom), dom, give_ope(dom, p1), n_cmp, p_cmp) ;
		
		else 
			eq[neq] = new Eq_one_side(espace.get_domain(dom), dom, 
						new Ope_sub(this, give_ope(dom, p1), give_ope(dom, p2)), n_cmp, p_cmp) ;
		neq ++ ;
	}
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_one_side (int dom, const char* nom, const List_comp& list) {
	add_eq_one_side (dom, nom, list.get_ncomp(), list.get_pcomp()) ;
}

 void System_of_eqs::add_eq_mode (int dom, int bound, const char* nom, const Index& pos_cf, double value) {

	char auxi[LMAX] ;
	trim_spaces(auxi, nom) ;

	eq_int[neq_int] = new Eq_int (1) ;
	eq_int[neq_int]->set_part(0, new Ope_mode(this, bound, pos_cf, value, give_ope(dom, auxi, bound))) ;

	neq_int ++ ;
	nbr_conditions = -1 ;
}

 void System_of_eqs::add_eq_val_mode (int dom, const char* nom, const Index& pos_cf, double value) {

	char auxi[LMAX] ;
	trim_spaces(auxi, nom) ;

	eq_int[neq_int] = new Eq_int (1) ;
	eq_int[neq_int]->set_part(0, new Ope_val_mode(this, pos_cf, value, give_ope(dom, auxi))) ;

	neq_int ++ ;
	nbr_conditions = -1 ;
}

 void System_of_eqs::add_eq_val (int dom, const char* nom, const Index& pos) {

	char auxi[LMAX] ;
	trim_spaces(auxi, nom) ;

	eq_int[neq_int] = new Eq_int (1) ;
	eq_int[neq_int]->set_part(0, new Ope_val(this, pos, give_ope(dom, auxi))) ;

	neq_int ++ ;
	nbr_conditions = -1 ;
}
void System_of_eqs::add_eq_point (int dom, const char* nom, const Point& num) {

	char auxi[LMAX] ;
	trim_spaces(auxi, nom) ;

	eq_int[neq_int] = new Eq_int (1) ;
	eq_int[neq_int]->set_part(0, new Ope_point(this, num, give_ope(dom, auxi))) ;

	neq_int ++ ;
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_order (int dom, const Array<int>& order, const char* nom, int n_cmp, Array<int>** p_cmp) {
    // Is it written like =0 ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;
		
		// no lhs :
		if (indic)
			eq[neq] = new Eq_order_array(espace.get_domain(dom), dom, order, give_ope(dom, p1), n_cmp, p_cmp) ;
		
		else 
			eq[neq] = new Eq_order_array(espace.get_domain(dom), dom, order,
						new Ope_sub(this, give_ope(dom, p1), give_ope(dom, p2)), n_cmp, p_cmp) ;
		neq ++ ;
	}
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_order (int dom, const Array<int>& order, const char* nom, const List_comp& list) {
	add_eq_order (dom, order, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_bc (int dom, int bound, const Array<int>& order, const char* nom, int n_cmp, Array<int>** p_cmp) {
    // Is it written like =0 ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for boundary conditions" << endl ;
		abort() ;
	}
	else {
		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;
			
		// no lhs :
		if (indic)
			eq[neq] = new Eq_bc_order_array(espace.get_domain(dom), dom, bound, order, give_ope(dom, p1, bound), n_cmp, p_cmp) ;
		else
			eq[neq] = new Eq_bc_order_array(espace.get_domain(dom), dom, 
				bound, order, new Ope_sub(this, give_ope(dom, p1, bound), give_ope(dom, p2, bound)), n_cmp, p_cmp) ;

		neq ++ ;
	}	
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_bc (int dom, int bound, const Array<int>& order, const char* nom, const List_comp& list) {
	add_eq_bc (dom, bound, order, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_matching (int dom, int bound, const Array<int>& order, const char* nom, int n_cmp, Array<int>** p_cmp) {
      int other_dom ;
	int other_bound ;
	espace.get_domain(dom)->find_other_dom (dom, bound, other_dom, other_bound) ;
	assert (other_dom>=dom_min) ;
	assert (other_dom<=dom_max) ;

	// Is it written with =  ?
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;

	if (!indic) {
		char auxi[LMAX] ;
		trim_spaces(auxi, nom) ;
		// Version without =
		eq[neq] = new Eq_matching_order_array(espace.get_domain(dom), dom, bound, other_dom, other_bound, order, 
			give_ope(dom, auxi, bound), give_ope(other_dom, auxi, other_bound), n_cmp, p_cmp) ;
		neq++ ;
	}
	else {
		// Version with = 
		eq[neq] = new Eq_matching_order_array(espace.get_domain(dom), dom, bound, other_dom, other_bound, order, 
			give_ope(dom, p1, bound), give_ope(other_dom, p2, other_bound), n_cmp, p_cmp) ;
		neq++ ;
	}	
	nbr_conditions = -1 ;
}

void System_of_eqs::add_eq_matching (int dom, int bound, const Array<int>& order, const char* nom, const List_comp& list) {
	add_eq_matching (dom, bound, order, nom, list.get_ncomp(), list.get_pcomp()) ;
}

void System_of_eqs::add_eq_first_integral (int dom_min, int dom_max, const char* integ_part, const char* cst_part) {

	eq[neq] = new Eq_first_integral(this, espace.get_domain(dom_min), dom_min, dom_max, integ_part, cst_part) ;
	neq++ ;

	nbr_conditions = -1 ;
}

}
