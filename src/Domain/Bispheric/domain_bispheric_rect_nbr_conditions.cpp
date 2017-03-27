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
int Domain_bispheric_rect::nbr_conditions_val_domain (const Val_domain& so, int order) const {
	
	int forgot_chi = 0;
	switch (order) {
	  case 0 :
	      forgot_chi = 0 ;
	      break ;
	  case 1 :
	      forgot_chi = 1 ;
	      break ;
	  case 2 :
	      forgot_chi = 1 ;
	      break ;
	  default:
	      cerr << "Unknown order in Domain_bispheric_rect::nbr_conditions_val_domain" << endl ;
	      abort() ;
	}
	
	int res = 0 ;
	int basep = (*so.get_base().bases_1d[2]) (0) ;
	
	// Loop on phi :
	for (int k=0 ; k<nbr_coefs(2) ; k++)
		// Loop on chi ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
			bool true_other = true ;
			
			switch (basep) {
				case COS :
				// Last odd ones
				if ((k%2==1) && (j==nbr_coefs(1)-1-forgot_chi))
					true_other = false ;
				// Regularity for even ones :
				if ((k!=0) && (k%2==0) && (j==0))
					true_other = false ;
				if (j==nbr_coefs(1)-forgot_chi)
					true_other = false ;
				break ;
			case SIN :
				// sin(0)
				if ((k==0) || (k==nbr_coefs(2)-1))
					true_other = false ;
				// Last odd ones :
				if ((k%2==1) && (j==nbr_coefs(1)-1-forgot_chi))
					true_other = false ;
				// Regularity for even ones :
				if ((k%2==0) && (j==0))
					true_other = false ;
				if (j==nbr_coefs(1)-forgot_chi)
					true_other = false ;
				break ;
			default :
				cerr << "Unknwon phi basis in Domain_bispheric_rect:nbr_conditions" << endl ;
				abort() ;
			}
			
			if (true_other)
				res += nbr_coefs(0)-order ;
	}
	return res ;
}

Array<int> Domain_bispheric_rect::nbr_conditions (const Tensor& tt, int dom, int order, int n_cmp, Array<int>** p_cmp) const {

	int size = (n_cmp==-1) ? tt.get_n_comp() : n_cmp ;
	Array<int> res (size) ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			res.set(0) = nbr_conditions_val_domain (tt()(dom), order) ;
			break ;
		case 1 : {
			bool found = false ;
			// Cartesian basis
			if (tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) {
				if (n_cmp==-1) {
					res.set(0) = nbr_conditions_val_domain (tt(1)(dom), order) ;
					res.set(1) = nbr_conditions_val_domain (tt(2)(dom), order) ;
					res.set(2) = nbr_conditions_val_domain (tt(3)(dom), order) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if ((*p_cmp[i])(0)==1)
						res.set(i) = nbr_conditions_val_domain (tt(1)(dom), order) ;
					if ((*p_cmp[i])(0)==2)
						res.set(i) = nbr_conditions_val_domain (tt(2)(dom), order) ;
					if ((*p_cmp[i])(0)==3)
						res.set(i) = nbr_conditions_val_domain (tt(3)(dom), order) ;
				}
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of vector Domain_bispheric_rect::nbr_conditions" << endl ;
				abort() ;
			}
		}
		      break ;
		case 2 : {
			bool found = false ;
			// Cartesian basis and symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==6)) {
				if (n_cmp==-1) {
					res.set(0) = nbr_conditions_val_domain (tt(1,1)(dom), order) ;
					res.set(1) = nbr_conditions_val_domain (tt(1,2)(dom), order) ;
					res.set(2) = nbr_conditions_val_domain (tt(1,3)(dom), order) ;
					res.set(3) = nbr_conditions_val_domain (tt(2,2)(dom), order) ;
					res.set(4) = nbr_conditions_val_domain (tt(2,3)(dom), order) ;
					res.set(5) = nbr_conditions_val_domain (tt(3,3)(dom), order) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==1))
						res.set(i) = nbr_conditions_val_domain (tt(1, 1)(dom), order) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==2))
						res.set(i) = nbr_conditions_val_domain (tt(1, 2)(dom), order) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(1, 3)(dom), order) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==2))
						res.set(i) = nbr_conditions_val_domain (tt(2, 2)(dom), order) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(2, 3)(dom), order) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(3, 3)(dom), order) ;
				}
				found = true ;
			}
			// Cartesian basis and not symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==9)) {
				if (n_cmp==-1) {
					res.set(0) = nbr_conditions_val_domain (tt(1,1)(dom), order) ;
					res.set(1) = nbr_conditions_val_domain (tt(1,2)(dom), order) ;
					res.set(2) = nbr_conditions_val_domain (tt(1,3)(dom), order) ;
					res.set(3) = nbr_conditions_val_domain (tt(2,1)(dom), order) ;
					res.set(4) = nbr_conditions_val_domain (tt(2,2)(dom), order) ;
					res.set(5) = nbr_conditions_val_domain (tt(2,3)(dom), order) ;
					res.set(6) = nbr_conditions_val_domain (tt(3,1)(dom), order) ;
					res.set(7) = nbr_conditions_val_domain (tt(3,2)(dom), order) ;
					res.set(8) = nbr_conditions_val_domain (tt(3,3)(dom),  order) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==1))
						res.set(i) = nbr_conditions_val_domain (tt(1, 1)(dom), order) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==2))
						res.set(i) = nbr_conditions_val_domain (tt(1, 2)(dom), order) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(1, 3)(dom), order) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==1))
						res.set(i) = nbr_conditions_val_domain (tt(2, 1)(dom), order) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==2))
						res.set(i) = nbr_conditions_val_domain (tt(2, 2)(dom), order) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(2, 3)(dom), order) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==1))
						res.set(i) = nbr_conditions_val_domain (tt(3, 1)(dom), order) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==2))
						res.set(i) = nbr_conditions_val_domain (tt(3, 2)(dom), order) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==3))
						res.set(i) = nbr_conditions_val_domain (tt(3, 3)(dom), order) ;
				}
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of 2-tensor Domain_bispheric_rect::nbr_conditions" << endl ;
				abort() ;
			}
		}
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_bispheric_rect::nbr_conditions" << endl ;
			break ;
	}
	return res ;
}}
