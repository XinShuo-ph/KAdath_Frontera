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

#include "headcpp.hpp"

#include "bispheric.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"
namespace Kadath {
void Space_bispheric::add_bc_sphere_one (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
  
	if (ndom_minus==1) {
	  // BC in domain 2 and 3
	  sys.add_eq_bc (ndom_minus+ndom_plus, INNER_BC, name, nused, pused) ;
	  sys.add_eq_bc (ndom_minus+ndom_plus+1, INNER_BC, name, nused, pused) ;
	}
	else {
	  sys.add_eq_bc (2, INNER_BC, name, nused, pused) ;
	}
}

void Space_bispheric::add_bc_sphere_two (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
	if (ndom_plus==1) {
	  // BC in domain 5 and 6
	  sys.add_eq_bc (ndom_minus+ndom_plus+3, INNER_BC, name, nused, pused) ;
	  sys.add_eq_bc (ndom_minus+ndom_plus+4, INNER_BC, name, nused, pused) ;
	}
	else {
	    sys.add_eq_bc (ndom_minus+1, INNER_BC, name, nused, pused) ;
	}
}

void Space_bispheric::add_bc_outer_sphere (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
	// All domains :	
	for (int dd=ndom_minus+ndom_plus ; dd<ndom_minus+ndom_plus+5 ; dd++)
		sys.add_eq_bc(dd, OUTER_BC, name, nused, pused) ;
}

void Space_bispheric::add_bc_outer (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {	
	sys.add_eq_bc (nbr_domains-1, OUTER_BC, name, nused, pused) ;
}

void Space_bispheric::add_eq (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {	
	for (int dd=sys.get_dom_min() ; dd<=sys.get_dom_max() ; dd++)
		sys.add_eq_inside (dd, name, nused, pused) ;
}

void Space_bispheric::add_eq_full (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {	
	for (int dd=sys.get_dom_min() ; dd<=sys.get_dom_max() ; dd++)
		sys.add_eq_full (dd, name, nused, pused) ;
}

void Space_bispheric::add_eq_one_side (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {	
	for (int dd=sys.get_dom_min() ; dd<=sys.get_dom_max() ; dd++)
		sys.add_eq_one_side (dd, name, nused, pused) ;
}

void Space_bispheric::add_matching (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
	sys.add_eq_matching (ndom_minus+ndom_plus, CHI_ONE_BC, name, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus+1, ETA_PLUS_BC, name, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus+2, ETA_PLUS_BC, name, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus+3, CHI_ONE_BC, name, nused, pused) ;
}

void Space_bispheric::add_matching_one_side (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
	sys.add_eq_matching_one_side (ndom_minus+ndom_plus, CHI_ONE_BC, name, nused, pused) ;
	sys.add_eq_matching_one_side (ndom_minus+ndom_plus+1, ETA_PLUS_BC, name, nused, pused) ;
	sys.add_eq_matching_one_side (ndom_minus+ndom_plus+2, ETA_PLUS_BC, name, nused, pused) ;
	sys.add_eq_matching_one_side (ndom_minus+ndom_plus+3, CHI_ONE_BC, name, nused, pused) ;
}

void Space_bispheric::add_eq (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused)  {
  
	// First sphere if any:
	for (int i=2 ; i<ndom_minus ; i++) {
	  sys.add_eq_inside (i, eq, nused, pused) ;
	  sys.add_eq_matching (i, OUTER_BC, rac, nused, pused) ;
	  sys.add_eq_matching (i, OUTER_BC, rac_der, nused, pused) ;
	}
	
	if (ndom_minus>1) {
	  sys.add_eq_inside (ndom_minus, eq, nused, pused) ;
	  // Matching with bispheric
	  sys.add_eq_matching_import (ndom_minus, OUTER_BC, rac, nused, pused) ;
	  sys.add_eq_matching_import (ndom_minus+ndom_plus, INNER_BC, rac_der, nused, pused) ;  
	  sys.add_eq_matching_import (ndom_minus+ndom_plus+1, INNER_BC, rac_der, nused, pused) ;
	}
	
	// second sphere (if any)
	for (int i=ndom_minus+1 ; i<ndom_minus+ndom_plus-1 ; i++) {
	  sys.add_eq_inside (i, eq, nused, pused) ;
	  sys.add_eq_matching (i, OUTER_BC, rac, nused, pused) ;
	  sys.add_eq_matching (i, OUTER_BC, rac_der, nused, pused) ;
	}
	
	if (ndom_plus>1) {
	  sys.add_eq_inside (ndom_minus+ndom_plus-1, eq, nused, pused) ;
	  // Matching with bispheric
	  sys.add_eq_matching_import (ndom_minus+ndom_plus-1, OUTER_BC, rac, nused, pused) ;
	  sys.add_eq_matching_import (ndom_minus+ndom_plus+3, INNER_BC, rac_der, nused, pused) ;
	  sys.add_eq_matching_import (ndom_minus+ndom_plus+4, INNER_BC, rac_der, nused, pused) ;
	}
	
	// Chi first
	sys.add_eq_inside (ndom_minus+ndom_plus, eq, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus, CHI_ONE_BC, rac_der, nused, pused) ;

	// Rect :
	sys.add_eq_inside (ndom_minus+ndom_plus+1, eq, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus+1, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus+1, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Eta first
	sys.add_eq_inside (ndom_minus+ndom_plus+2, eq, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus+2, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus+2, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Rect 
	sys.add_eq_inside (ndom_minus+ndom_plus+3, eq, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus+3, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (ndom_minus+ndom_plus+3, CHI_ONE_BC, rac_der, nused, pused) ;

	// chi first :
	sys.add_eq_inside (ndom_minus+ndom_plus+4, eq, nused, pused) ;

	// Matching outer domain :
	if (ndom_minus+ ndom_plus + 5 <= sys.get_dom_max()) {
	for (int i=ndom_minus+ndom_plus ; i<=ndom_minus+ndom_plus+4 ; i++)
		sys.add_eq_matching_import (i, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (ndom_minus+ndom_plus+5, INNER_BC, rac_der, nused, pused) ;
	
	// Shells :
	for (int i=0 ; i<nshells ; i++) {
	    sys.add_eq_inside (ndom_minus+ndom_plus+5+i, eq, nused, pused) ;
	    sys.add_eq_matching (ndom_minus+ndom_plus+6+i, INNER_BC, rac, nused, pused) ;
	    sys.add_eq_matching (ndom_minus+ndom_plus+6+i, INNER_BC, rac_der, nused, pused) ;
	}

	 //Compactified domain
	sys.add_eq_inside (nbr_domains-1, eq, nused, pused) ;}
}

void Space_bispheric::add_eq (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, const List_comp& list) {
 add_eq (sys, eq, rac, rac_der, list.get_ncomp(), list.get_pcomp()) ;
}

void Space_bispheric::add_eq_int_inf (System_of_eqs& sys, const char* nom) {

	// Check the last domain is of the right type :
	const Domain_compact* pcomp = dynamic_cast <const Domain_compact*> (domains[nbr_domains-1]) ;
	if (pcomp==0x0) {
		cerr << "add_eq_int_inf requires a compactified domain" << endl ;
		abort() ;
	}
	int dom = nbr_domains-1 ;

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;

		// Construction of the equation
		sys.eq_int[sys.neq_int] = new Eq_int(1) ;

		// Affectation :
		// no lhs :
		if (indic)
			sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(dom, p1, OUTER_BC)) ;
		
		else 
			sys.eq_int[sys.neq_int]->set_part(0, new Ope_sub(&sys, sys.give_ope(dom, p1, OUTER_BC), sys.give_ope(dom, p2, OUTER_BC))) ;
		sys.neq_int ++ ;
	}
	sys.nbr_conditions = -1 ;
}

void Space_bispheric::add_eq_int_sphere_one (System_of_eqs& sys, const char* nom) {

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {

		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;

		if (ndom_minus==1) {
		  // Construction of the equation
		  sys.eq_int[sys.neq_int] = new Eq_int(2) ;

		  // Affectation :
		  // no lhs :
		  if (indic) {
			sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(ndom_minus+ndom_plus, p1, INNER_BC)) ;
			sys.eq_int[sys.neq_int]->set_part(1, sys.give_ope(ndom_minus+ndom_plus+1, p1, INNER_BC)) ;
			}
		  else {
			sys.eq_int[sys.neq_int]->set_part(0, 
				new Ope_sub(&sys, sys.give_ope(ndom_minus+ndom_plus, p1, INNER_BC), sys.give_ope(ndom_minus+ndom_plus, p2, INNER_BC))) ;
			sys.eq_int[sys.neq_int]->set_part(1, 
				new Ope_sub(&sys, sys.give_ope(ndom_minus+ndom_plus+1, p1, INNER_BC), sys.give_ope(ndom_minus+ndom_plus+1, p2, INNER_BC))) ;
		  }
		}
		else  {
		  // Construction of the equation
		  sys.eq_int[sys.neq_int] = new Eq_int(1) ;

		  // Affectation :
		  // no lhs :
		  if (indic) {
			sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(1, p1, INNER_BC)) ;
			}
		  else {
			sys.eq_int[sys.neq_int]->set_part(0, 
				new Ope_sub(&sys, sys.give_ope(1, p1, INNER_BC), sys.give_ope(1, p2, INNER_BC))) ;
		  }
		}
		sys.neq_int ++ ;
	}
	sys.nbr_conditions = -1 ;
}


void Space_bispheric::add_eq_int_sphere_two (System_of_eqs& sys, const char* nom) {

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;

		if (ndom_plus==1) {
		// Construction of the equation
		sys.eq_int[sys.neq_int] = new Eq_int(2) ;

		// Affectation :
		// no lhs :
		if (indic) {
			sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(ndom_minus+ndom_plus+3, p1, INNER_BC)) ;
			sys.eq_int[sys.neq_int]->set_part(1, sys.give_ope(ndom_minus+ndom_plus+4, p1, INNER_BC)) ;
			}
		else {
			sys.eq_int[sys.neq_int]->set_part(0, 
				new Ope_sub(&sys, sys.give_ope(ndom_minus+ndom_plus+3, p1, INNER_BC), sys.give_ope(ndom_minus+ndom_plus+3, p2, INNER_BC))) ;
			sys.eq_int[sys.neq_int]->set_part(1, 
				new Ope_sub(&sys, sys.give_ope(ndom_minus+ndom_plus+4, p1, INNER_BC), sys.give_ope(ndom_minus+ndom_plus+4, p2, INNER_BC))) ;
		}
		}
		else {
		  // Construction of the equation
		sys.eq_int[sys.neq_int] = new Eq_int(1) ;

		// Affectation :
		// no lhs :
		if (indic) {
			sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(ndom_minus+1, p1, INNER_BC)) ;
			}
		else {
			sys.eq_int[sys.neq_int]->set_part(0, 
				new Ope_sub(&sys, sys.give_ope(ndom_minus+1, p1, INNER_BC), sys.give_ope(ndom_minus+1, p2, INNER_BC))) ;
		}
		}
		  
		sys.neq_int ++ ;
	}
	sys.nbr_conditions = -1 ;
}


void Space_bispheric::add_eq_zero_mode_inf (System_of_eqs& sys, const char* name, int j, int k) {

	Index pos_cf (domains[nbr_domains-1]->get_nbr_coefs()) ;
	pos_cf.set(1) = j ;
	pos_cf.set(2) = k ;
	double value = 0. ;
	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	sys.add_eq_mode (nbr_domains-1, OUTER_BC, auxi, pos_cf, value) ;
}
}
