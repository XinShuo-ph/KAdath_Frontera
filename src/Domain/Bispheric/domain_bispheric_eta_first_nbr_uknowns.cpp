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
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
int Domain_bispheric_eta_first::nbr_unknowns_val_domain (const Val_domain& so) const {
	
	int res = 0 ;
	int basep = (*so.get_base().bases_1d[2]) (0) ;

	Index pos (nbr_coefs) ;
	do {
		bool indic = true ;
		switch (basep) {
		case COS :
			// Last odd ones
			if ((pos(2)%2==1) && (pos(0)==nbr_coefs(0)-1))
				indic = false ;
			// Regularity for even ones :
			if ((pos(2)!=0) && (pos(2)%2==0) && (pos(0)==0))
				indic = false ;
			break ;
		case SIN :
			// sin(0)
			if ((pos(2)==0) || (pos(2)==nbr_coefs(2)-1))
				indic = false ;
			// Last odd ones :
			if ((pos(2)%2==1) && (pos(0)==nbr_coefs(0)-1))
				indic = false ;
			// Regularity for even ones :
			if ((pos(2)%2==0) && (pos(0)==0))
				indic = false ;
			break ;
		default :
			cerr << "Unknown phi basis in Domain_bispheric_eta_first::nbr_unknwowns_val_domain" << endl ;
			abort() ;
		}
		if (indic)
			res ++ ;
	}
	while (pos.inc()) ;

	return res ;
}

int Domain_bispheric_eta_first::nbr_unknowns (const Tensor& tt, int dom) const {

	// Check right domain
	assert (tt.get_space().get_domain(dom)==this) ;

	int res = 0 ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			res += nbr_unknowns_val_domain (tt()(dom)) ;
			break ;
		case 1 : {
			bool found = false ;
			// Cartesian basis
			if (tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) {
				res += nbr_unknowns_val_domain (tt(1)(dom)) ;
				res += nbr_unknowns_val_domain (tt(2)(dom)) ;
				res += nbr_unknowns_val_domain (tt(3)(dom)) ;
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of vector Domain_bispheric_eta_first::nbr_unknowns" << endl ;
				abort() ;
			}
		}
		    break ;
		case 2 : {
			bool found = false ;
			// Cartesian basis and symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==6)) {
				res += nbr_unknowns_val_domain (tt(1,1)(dom)) ;
				res += nbr_unknowns_val_domain (tt(1,2)(dom)) ;
				res += nbr_unknowns_val_domain (tt(1,3)(dom)) ;
				res += nbr_unknowns_val_domain (tt(2,2)(dom)) ;
				res += nbr_unknowns_val_domain (tt(2,3)(dom)) ;
				res += nbr_unknowns_val_domain (tt(3,3)(dom)) ;
				found = true ;
			}
			// Cartesian basis and not symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==9)) {
				for (int i=1 ; i<=3 ; i++)
				  for (int j=1 ; j<=3 ; j++)
				    res += nbr_unknowns_val_domain (tt(i,j)(dom)) ;
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of 2-tensor Domain_bispheric_eta_first::nbr_unknowns" << endl ;
				abort() ;
			}
		}
		break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_bispheric_eta_first::nbr_unknowns" << endl ;
			break ;
	}
	return res ;
}
}
