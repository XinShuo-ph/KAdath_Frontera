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
#include "polar.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"
#include "ope_eq.hpp"
namespace Kadath {
void Space_polar::add_inner_bc (System_of_eqs& sys, const char* name, int nused, Array<int>** pused) const   {
	// BC in domain 5 and 6
	sys.add_eq_bc (sys.get_dom_min(), INNER_BC, name, nused, pused) ;
}

void Space_polar::add_outer_bc (System_of_eqs& sys, const char* name, int nused, Array<int>** pused) const  {
	sys.add_eq_bc (sys.get_dom_max(), OUTER_BC, name, nused, pused) ;
}

void Space_polar::add_eq (System_of_eqs& sys, const char* name, int nused, Array<int>** pused) const {
	for (int dd=sys.get_dom_min() ; dd<=sys.get_dom_max() ; dd++)
		sys.add_eq_inside (dd, name, nused, pused) ;
}

void Space_polar::add_eq_full (System_of_eqs& sys, const char* name, int nused, Array<int>** pused) const {
	for (int dd=sys.get_dom_min() ; dd<=sys.get_dom_max() ; dd++)
		sys.add_eq_full (dd, name, nused, pused) ;
}

void Space_polar::add_matching (System_of_eqs& sys, const char* name, int nused, Array<int>** pused) const {
	for (int dd=sys.get_dom_min() ; dd<sys.get_dom_max() ; dd++)
		sys.add_eq_matching (dd, OUTER_BC, name, nused, pused) ;
}

void Space_polar::add_eq (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused) const  {
	for (int dd=sys.get_dom_min() ; dd<sys.get_dom_max(); dd++) {
		sys.add_eq_inside (dd, eq, nused, pused) ;
		sys.add_eq_matching (dd, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (dd, OUTER_BC, rac_der, nused, pused) ;
	}
	sys.add_eq_inside (sys.get_dom_max(), eq, nused, pused) ;
}

void Space_polar::add_eq_int_volume (System_of_eqs& sys, const char* nom) {
  sys.eq_int_list.push_back(std::make_tuple(nom, nbr_domains-1, -1));

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
		sys.eq_int[sys.neq_int] = new Eq_int(nbr_domains+1) ;

		// Affectation of the intregrale parts
		for (int d=0 ; d<nbr_domains ; d++)
		  sys.eq_int[sys.neq_int]->set_part(d, sys.give_ope(d, p1)) ;
		// Affectation of the second member (constant value)
		sys.eq_int[sys.neq_int]->set_part(nbr_domains, new Ope_minus(&sys, sys.give_ope(0, p2))) ;
		sys.neq_int ++ ;
	}
	sys.nbr_conditions = -1 ;
}


void Space_polar::add_eq_int_inf (System_of_eqs& sys, const char* nom) {
  sys.eq_int_list.push_back(std::make_tuple(nom, nbr_domains-1, OUTER_BC));

	// Check the last domain is of the right type :
	const Domain_polar_compact* pcomp = dynamic_cast <const Domain_polar_compact*> (domains[nbr_domains-1]) ;
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

void Space_polar::add_eq_int_inner (System_of_eqs& sys, const char* nom) {
  sys.eq_int_list.push_back(std::make_tuple(nom, 1, OUTER_BC));

	// Check the last domain is of the right type :
	const Domain_polar_shell* pshell = dynamic_cast <const Domain_polar_shell*> (domains[1]) ;
	if (pshell==0x0) {
		cerr << "add_eq_int_inner requires that the domain 1 is a shell" << endl ;
		abort() ;
	}
	int dom = 1 ;

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


void Space_polar::add_eq_mode (System_of_eqs& sys, const char* name, int domtarget, int itarget, int jtarget) {

	Index pos_cf (domains[domtarget]->get_nbr_coefs()) ;
	pos_cf.set(0) = itarget ;
	pos_cf.set(1) = jtarget ;
	double value = 1. ;
	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	sys.add_eq_val_mode (domtarget, auxi, pos_cf, value) ;
}

void Space_polar::add_eq_ori (System_of_eqs& sys, const char* name) {

	Index pos (domains[0]->get_nbr_points()) ;
	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	sys.add_eq_val (0, auxi, pos) ;
}

void Space_polar::add_eq_point (System_of_eqs& sys, const Point& MM, const char* name) {

	// Get the domain and num coordinates :
	int ld = -1;
	bool* inside = new bool[nbr_domains] ;
	for (int l=0 ; l<nbr_domains; l++)
	    inside[l] = get_domain(l)->is_in(MM) ;
	for (int l=0 ; l<nbr_domains ; l++)
	      if ((ld==-1) && (inside[l]))
		  ld = l ;
	Point num (get_domain(ld)->absol_to_num(MM)) ;

	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	sys.add_eq_point (ld, auxi, num) ;
}
}
