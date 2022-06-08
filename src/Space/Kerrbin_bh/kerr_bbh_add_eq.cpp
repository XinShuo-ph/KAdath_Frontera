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
#include "kerrbin_bh.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"
namespace Kadath {
void Space_Kerr_bbh::add_bc_sphere_one (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
	  sys.add_eq_bc (BH1+2, INNER_BC, name, nused, pused) ;
}

void Space_Kerr_bbh::add_bc_sphere_two (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
	sys.add_eq_bc (BH2+2, INNER_BC, name, nused, pused) ;
}

void Space_Kerr_bbh::add_bc_outer (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
	sys.add_eq_bc (OUTER+5, OUTER_BC, name, nused, pused) ;
}

void Space_Kerr_bbh::add_eq (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused)  {
	// First bh
	sys.add_eq_inside (BH1+2, eq, nused, pused) ;
  for(int i = 0; i < n_shells1; ++i) {
  	sys.add_eq_matching (BH1+2+i, OUTER_BC, rac, nused, pused) ;
  	sys.add_eq_matching (BH1+2+i, OUTER_BC, rac_der, nused, pused) ;
  	sys.add_eq_inside (BH1+3+i, eq, nused, pused) ;
  }

	// Matching with bispheric :
	sys.add_eq_matching_import (BH1 + n_shells1 + 2, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (OUTER, INNER_BC, rac_der, nused, pused) ;
	sys.add_eq_matching_import (OUTER+1, INNER_BC, rac_der, nused, pused) ;

	// Second BH :
	sys.add_eq_inside (BH2+2, eq, nused, pused) ;
  for(int i = 0; i < n_shells2; ++i) {
  	sys.add_eq_matching (BH2+2+i, OUTER_BC, rac, nused, pused) ;
  	sys.add_eq_matching (BH2+2+i, OUTER_BC, rac_der, nused, pused) ;
  	sys.add_eq_inside (BH2+3+i, eq, nused, pused) ;
  }

	// Matching with bispheric :
	sys.add_eq_matching_import (BH2 + n_shells2 + 2, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (OUTER+3, INNER_BC, rac_der, nused, pused) ;
	sys.add_eq_matching_import (OUTER+4, INNER_BC, rac_der, nused, pused) ;

	// Chi first
	sys.add_eq_inside (OUTER, eq, nused, pused) ;
	sys.add_eq_matching (OUTER, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER, CHI_ONE_BC, rac_der, nused, pused) ;

	// Rect :
	sys.add_eq_inside (OUTER+1, eq, nused, pused) ;
	sys.add_eq_matching (OUTER+1, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER+1, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Eta first
	sys.add_eq_inside (OUTER+2, eq, nused, pused) ;
	sys.add_eq_matching (OUTER+2, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER+2, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Rect
	sys.add_eq_inside (OUTER+3, eq, nused, pused) ;
	sys.add_eq_matching (OUTER+3, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER+3, CHI_ONE_BC, rac_der, nused, pused) ;


	// chi first :
	sys.add_eq_inside (OUTER+4, eq, nused, pused) ;

	// Matching outer domain :
	for (int d=OUTER ; d<=OUTER+4 ; d++)
		sys.add_eq_matching_import (d, OUTER_BC, rac, nused, pused) ;
  
  //Matching for first shell or compactified domain
	sys.add_eq_matching_import (OUTER+5, INNER_BC, rac_der, nused, pused) ;
  
  // Optional spherical shells between bi-spherical and compactified domain
  for (int d=0 ; d<n_shells_outer ; d++) {
      sys.add_eq_inside (OUTER+5+d, eq, nused, pused) ;
      sys.add_eq_matching (OUTER+5+d, OUTER_BC, rac, nused, pused) ;
      sys.add_eq_matching (OUTER+5+d, OUTER_BC, rac_der, nused, pused) ;
  }

  //Compactified domain
  sys.add_eq_inside (OUTER+5+n_shells_outer, eq, nused, pused) ;

}

void Space_Kerr_bbh::add_eq_int_inf (System_of_eqs& sys, const char* nom) {

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
    sys.eq_int_list.push_back(std::make_tuple(nom, nbr_domains-1, OUTER_BC));
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


void Space_Kerr_bbh::add_eq_int_sphere_one (System_of_eqs& sys, const char* nom) {

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
    sys.eq_int_list.push_back(std::make_tuple(nom,BH1+2,INNER_BC));

		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;

	  // Construction of the equation
	  sys.eq_int[sys.neq_int] = new Eq_int(1) ;

		  // Affectation :
		  // no lhs :
		  if (indic) {
			sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(BH1+2, p1, INNER_BC)) ;
			}
		  else {
			sys.eq_int[sys.neq_int]->set_part(0,
				new Ope_sub(&sys, sys.give_ope(BH1+2, p1, INNER_BC), sys.give_ope(BH1+2, p2, INNER_BC))) ;
		  }
		}
		sys.neq_int ++ ;
	sys.nbr_conditions = -1 ;
}

void Space_Kerr_bbh::add_eq_int_sphere_two (System_of_eqs& sys, const char* nom) {

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
    sys.eq_int_list.push_back(std::make_tuple(nom,BH2+2,INNER_BC));

		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;

	  // Construction of the equation
	  sys.eq_int[sys.neq_int] = new Eq_int(1) ;

		  // Affectation :
		  // no lhs :
		  if (indic) {
			sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(BH2+2, p1, INNER_BC)) ;
			}
		  else {
			sys.eq_int[sys.neq_int]->set_part(0,
				new Ope_sub(&sys, sys.give_ope(BH2+2, p1, INNER_BC), sys.give_ope(BH2+2, p2, INNER_BC))) ;
		  }
		}
		sys.neq_int ++ ;
	sys.nbr_conditions = -1 ;
}

void Space_Kerr_bbh::add_eq_zero_mode_inf (System_of_eqs& sys, const char* name, int j, int k) {
	Index pos_cf (domains[nbr_domains-1]->get_nbr_coefs()) ;
	pos_cf.set(1) = j ;
	pos_cf.set(2) = k ;
	double value = 0. ;
	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	sys.add_eq_mode (nbr_domains-1, OUTER_BC, auxi, pos_cf, value) ;
}

}
