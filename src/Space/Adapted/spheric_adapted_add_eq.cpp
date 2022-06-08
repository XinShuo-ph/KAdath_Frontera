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
#include "adapted.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"
namespace Kadath {
void Space_spheric_adapted::add_eq_int (System_of_eqs& sys, const int dom, const int bc, const char* eq) {

  sys.eq_int_list.push_back(std::make_tuple(eq,dom,bc));

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(eq, p1, p2, '=') ;
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
		  if (indic) {
	  		sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(dom, p1, bc)) ;
			}
		  else {
		  	sys.eq_int[sys.neq_int]->set_part(0,
			  	new Ope_sub(&sys, sys.give_ope(dom, p1, bc), sys.give_ope(dom, p2, bc))) ;
		  }

		sys.neq_int ++ ;
	  sys.nbr_conditions = -1 ;
	}
}
void Space_spheric_adapted::add_eq_int_inf (System_of_eqs& sys, const char* nom) {

	// Check the last domain is of the right type :
	const Domain_compact* pcomp = dynamic_cast <const Domain_compact*> (domains[nbr_domains-1]) ;
	if (pcomp==0x0) {
		cerr << "add_eq_int_inf requires a compactified domain" << endl ;
		abort() ;
	}
	int dom = nbr_domains-1 ;
  sys.eq_int_list.push_back(std::make_tuple(nom,dom,OUTER_BC));

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Construction of the equation
		sys.eq_int[sys.neq_int] = new Eq_int(1) ;

		// Affectation :
		sys.eq_int[sys.neq_int]->set_part(0, new Ope_sub(&sys, sys.give_ope(dom, p1, OUTER_BC), sys.give_ope(dom, p2, OUTER_BC))) ;
		sys.neq_int ++ ;
	}
	sys.nbr_conditions = -1 ;
}

void Space_spheric_adapted::add_eq_int_volume (System_of_eqs& sys, int nz, const char* nom) {

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
    sys.eq_int_list.push_back(std::make_tuple(nom, 0, -1));
		// Construction of the equation
		sys.eq_int[sys.neq_int] = new Eq_int(nz+1) ;

		// Affectation of the intregrale parts
		for (int d=0 ; d<nz ; d++)
		  sys.eq_int[sys.neq_int]->set_part(d, sys.give_ope(d, p1)) ;
		// Affectation of the second member (constant value)
		sys.eq_int[sys.neq_int]->set_part(nz, new Ope_minus(&sys, sys.give_ope(0, p2))) ;
		sys.neq_int ++ ;
	}
	sys.nbr_conditions = -1 ;
}

void Space_spheric_adapted::add_eq (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused)  {
	for (int dd=sys.get_dom_min() ; dd<sys.get_dom_max(); dd++) {
		sys.add_eq_inside (dd, eq, nused, pused) ;
		sys.add_eq_matching (dd, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (dd, OUTER_BC, rac_der, nused, pused) ;
	}
	sys.add_eq_inside (sys.get_dom_max(), eq, nused, pused) ;
}}
