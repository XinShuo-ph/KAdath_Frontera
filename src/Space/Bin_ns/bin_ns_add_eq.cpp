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
#include "bin_ns.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"

namespace Kadath {
void Space_bin_ns::add_eq (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused)  {
  
	// First NS 
	sys.add_eq_inside (0, eq, nused, pused) ;
	sys.add_eq_matching (0, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (0, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (1, eq, nused, pused) ;
	sys.add_eq_matching (1, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (1, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (2, eq, nused, pused) ;
	
	
	// Matching with bispheric :
	sys.add_eq_matching_import (2, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (6, INNER_BC, rac_der, nused, pused) ;  
	sys.add_eq_matching_import (7, INNER_BC, rac_der, nused, pused) ;
	
	// Second NS :
	sys.add_eq_inside (3, eq, nused, pused) ;
	sys.add_eq_matching (3, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (3, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (4, eq, nused, pused) ;
	sys.add_eq_matching (4, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (4, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (5, eq, nused, pused) ;
	
	// Matching with bispheric :
	sys.add_eq_matching_import (5, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (9, INNER_BC, rac_der, nused, pused) ;  
	sys.add_eq_matching_import (10, INNER_BC, rac_der, nused, pused) ;
	
	// Chi first
	sys.add_eq_inside (6, eq, nused, pused) ;
	sys.add_eq_matching (6, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (6, CHI_ONE_BC, rac_der, nused, pused) ;

	// Rect :
	sys.add_eq_inside (7, eq, nused, pused) ;
	sys.add_eq_matching (7, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (7, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Eta first
	sys.add_eq_inside (8, eq, nused, pused) ;
	sys.add_eq_matching (8, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (8, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Rect 
	sys.add_eq_inside (9, eq, nused, pused) ;
	sys.add_eq_matching (9, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (9, CHI_ONE_BC, rac_der, nused, pused) ;

	
	// chi first :
	sys.add_eq_inside (10, eq, nused, pused) ;

	// Matching outer domain :
	for (int d=6 ; d<=10 ; d++)
		sys.add_eq_matching_import (d, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (11, INNER_BC, rac_der, nused, pused) ;
		
	 //Compactified domain
	sys.add_eq_inside (11, eq, nused, pused) ;
	
}

void Space_bin_ns::add_eq_nozec (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused)  {
  
	// First NS 
	sys.add_eq_inside (0, eq, nused, pused) ;
	sys.add_eq_matching (0, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (0, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (1, eq, nused, pused) ;
	sys.add_eq_matching (1, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (1, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (2, eq, nused, pused) ;
	
	
	// Matching with bispheric :
	sys.add_eq_matching_import (2, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (6, INNER_BC, rac_der, nused, pused) ;  
	sys.add_eq_matching_import (7, INNER_BC, rac_der, nused, pused) ;
	
	// Second NS :
	sys.add_eq_inside (3, eq, nused, pused) ;
	sys.add_eq_matching (3, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (3, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (4, eq, nused, pused) ;
	sys.add_eq_matching (4, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (4, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (5, eq, nused, pused) ;
	
	// Matching with bispheric :
	sys.add_eq_matching_import (5, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (9, INNER_BC, rac_der, nused, pused) ;  
	sys.add_eq_matching_import (10, INNER_BC, rac_der, nused, pused) ;
	
	// Chi first
	sys.add_eq_inside (6, eq, nused, pused) ;
	sys.add_eq_matching (6, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (6, CHI_ONE_BC, rac_der, nused, pused) ;

	// Rect :
	sys.add_eq_inside (7, eq, nused, pused) ;
	sys.add_eq_matching (7, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (7, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Eta first
	sys.add_eq_inside (8, eq, nused, pused) ;
	sys.add_eq_matching (8, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (8, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Rect 
	sys.add_eq_inside (9, eq, nused, pused) ;
	sys.add_eq_matching (9, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (9, CHI_ONE_BC, rac_der, nused, pused) ;

	
	// chi first :
	sys.add_eq_inside (10, eq, nused, pused) ;
	
	if (nshells !=0) {
	// Matching outer domain :
	for (int d=6 ; d<=10 ; d++)
		sys.add_eq_matching_import (d, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (11, INNER_BC, rac_der, nused, pused) ;

	for (int d=0 ; d<nshells ; d++) {
		sys.add_eq_inside (11+d, eq, nused, pused) ;
		if (d!=nshells-1) {
			sys.add_eq_matching (11+d, OUTER_BC, rac, nused, pused) ;
			sys.add_eq_matching (11+d, OUTER_BC, rac_der, nused, pused) ;
		}
	}
	}
}

void Space_bin_ns::add_eq_int_inf (System_of_eqs& sys, const char* nom) {

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

void Space_bin_ns::add_eq_ori_one (System_of_eqs& sys, const char* name) {

	Index pos (domains[0]->get_nbr_points()) ;
	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	sys.add_eq_val (0, auxi, pos) ;
}

void Space_bin_ns::add_eq_ori_two (System_of_eqs& sys, const char* name) {

	Index pos (domains[3]->get_nbr_points()) ;
	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	sys.add_eq_val (3, auxi, pos) ;
}}
