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
#include "point.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
int Domain_shell_outer_adapted::nbr_unknowns_val_domain (const Val_domain& so, int mlim) const {
	
	int res = 0 ;
	int kmin = 2*mlim + 2 ;

	Index pos (nbr_coefs) ;
	do {
		bool indic = true ;
		// True coef in phi ?
		if ((pos(2)==1) || (pos(2)==nbr_coefs(2)-1))
			indic = false ;
		// Get base in theta :
		int baset = (*so.get_base().bases_1d[1]) (pos(2)) ;
		switch (baset) {
					case COS_EVEN:
						if ((pos(1)==0) && (pos(2)>=kmin))
							indic = false ;
						break ;
					case COS_ODD:
						if ((pos(1)==nbr_coefs(1)-1) || ((pos(1)==0) && (pos(2)>=kmin)))
							indic = false ;
						break ;
					case SIN_EVEN:
						if (((pos(1)==1)&&(pos(2)>=kmin+2)) || (pos(1)==0) || (pos(1)==nbr_coefs(1)-1)) 
							indic = false  ;
						break ;
					case SIN_ODD:
						if (((pos(1)==0)&&(pos(2)>=kmin+2)) || (pos(1)==nbr_coefs(1)-1))
							indic = false ;
						break ;
					default:
						cerr << "Unknow theta basis in Domain_shell_outer_adapted::nbr_unknowns_val_domain" << endl ;
						abort() ;
		}
		if (indic)
			res ++ ;
	}
	while (pos.inc()) ;

	return res ;
}

int Domain_shell_outer_adapted::nbr_unknowns (const Tensor& tt, int dom) const {

	// Check right domain
	assert (tt.get_space().get_domain(dom)==this) ;

	int res = 0 ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			res += nbr_unknowns_val_domain (tt()(dom), 0) ;
			break ;
		case 1 : {
			bool found = false ;
			// Cartesian basis
			if (tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) {
				res += nbr_unknowns_val_domain (tt(1)(dom), 0) ;
				res += nbr_unknowns_val_domain (tt(2)(dom), 0) ;
				res += nbr_unknowns_val_domain (tt(3)(dom), 0) ;
				found = true ;
			}
			// Spherical coordinates
			if (tt.get_basis().get_basis(dom)==SPHERICAL_BASIS) {
				res += nbr_unknowns_val_domain (tt(1)(dom), 0) ;
				res += nbr_unknowns_val_domain (tt(2)(dom), 1) ;
				res += nbr_unknowns_val_domain (tt(3)(dom), 1) ;
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of vector Domain_shell_outer_adapted::nbr_unknowns" << endl ;
				abort() ;
			}
		}
		    break ;
		case 2 : {
			bool found = false ;
			// Cartesian basis and symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==6)) {
				res += nbr_unknowns_val_domain (tt(1,1)(dom), 0) ;
				res += nbr_unknowns_val_domain (tt(1,2)(dom), 0) ;
				res += nbr_unknowns_val_domain (tt(1,3)(dom), 0) ;
				res += nbr_unknowns_val_domain (tt(2,2)(dom), 0) ;
				res += nbr_unknowns_val_domain (tt(2,3)(dom), 0) ;
				res += nbr_unknowns_val_domain (tt(3,3)(dom), 0) ;
				found = true ;
			}
			// Cartesian basis and not symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==9)) {
				for (int i=1 ; i<=3 ; i++)
				  for (int j=1 ; j<=3 ; j++)
				    res += nbr_unknowns_val_domain (tt(i,j)(dom), 0) ;
				found = true ;
			}
			// Spherical coordinates and symetric
			if ((tt.get_basis().get_basis(dom)==SPHERICAL_BASIS) && (tt.get_n_comp()==6)) {
				res += nbr_unknowns_val_domain (tt(1,1)(dom), 0) ;
				res += nbr_unknowns_val_domain (tt(1,2)(dom), 1) ;
				res += nbr_unknowns_val_domain (tt(1,3)(dom), 1) ;
				res += nbr_unknowns_val_domain (tt(2,2)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(2,3)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(3,3)(dom), 2) ;
				found = true ;
			}
			// Spherical coordinates and not symetric
			if ((tt.get_basis().get_basis(dom)==SPHERICAL_BASIS) && (tt.get_n_comp()==9)) {
				res += nbr_unknowns_val_domain (tt(1,1)(dom), 0) ;
				res += nbr_unknowns_val_domain (tt(1,2)(dom), 1) ;
				res += nbr_unknowns_val_domain (tt(1,3)(dom), 1) ;
				res += nbr_unknowns_val_domain (tt(2,1)(dom), 1) ;
				res += nbr_unknowns_val_domain (tt(2,2)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(2,3)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(3,1)(dom), 1) ;
				res += nbr_unknowns_val_domain (tt(3,2)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(3,3)(dom), 2) ;
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of 2-tensor Domain_shell_outer_adapted::nbr_unknowns" << endl ;
				abort() ;
			}
		}
		break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_shell_outer_adapted::nbr_unknowns" << endl ;
			break ;
	}
	return res ;
}
}

