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
#include "polar_periodic.hpp"
#include "point.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
void Domain_polar_periodic_nucleus::export_tau_val_domain_boundary (const Val_domain& so, int bound, Array<double>& sec, int& pos_sec, int ncond) const {
	
	if (so.check_if_zero())
		pos_sec += ncond ;
	else {
	so.coef() ;
	Index pos_cf (nbr_coefs) ;
	
	int basetime = (*so.get_base().bases_1d[2]) (0) ;
	int mink, maxk ;
	switch (basetime) {
		case COS :
			mink=0 ;
			maxk=nbr_coefs(2) ;
			break ;
		case SIN :
			mink=1 ;
			maxk=nbr_coefs(2)-1 ;
			break ;
		default :
			cerr << "Unknown time basis in Domain_polar_periodic_nucleus_export_tau_val_domain_boundary" << endl ;
			abort() ;
	}


	for (int k=mink ; k<maxk ; k++) {
		pos_cf.set(2) = k ;

	// Loop on theta
	int baset = (*so.get_base().bases_1d[1]) (k) ;
	int minj, maxj ;
	switch (baset) {
		case COS_EVEN :
			minj=0 ;
			maxj=nbr_coefs(1) ;
			break ;
		case COS_ODD :
			minj=0 ;
			maxj=nbr_coefs(1)-1 ;
			break ;
		case SIN_EVEN :
			minj=1 ;
			maxj=nbr_coefs(1)-1 ;
			break ;
		case SIN_ODD :
			minj=0 ;
			maxj=nbr_coefs(1)-1 ;
			break ;
		
		default :
			cerr << "Unknown theta basis in Domain_polar_periodic_nucleus_export_tau_val_domain" << endl ;
			abort() ;
	}


	for (int j=minj ; j<maxj ; j++) {
		pos_cf.set(1) = j ;
		sec.set(pos_sec) = val_boundary(bound, so, pos_cf) ;
		pos_sec ++ ;
		} // end loop j 
	} // end loop k
	} // end null case
}

void Domain_polar_periodic_nucleus::export_tau_boundary (const Tensor& tt, int dom, int bound, Array<double>& res, int& pos_res, const Array<int>& ncond,
										int n_cmp, Array<int>** p_cmp) const {

	// Check boundary
	if (bound!=OUTER_BC) {
		cerr << "Unknown boundary in Domain_polar_periodic_nucleus::export_tau_boundary" << endl ;
		abort() ;
	}

	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			export_tau_val_domain_boundary (tt()(dom), bound, res, pos_res, ncond(0)) ;
			break ;
		case 1 :
			if (n_cmp==-1) {
					export_tau_val_domain_boundary (tt(1)(dom), bound, res, pos_res, ncond(0)) ;
					export_tau_val_domain_boundary (tt(2)(dom), bound,  res, pos_res, ncond(1)) ;
					export_tau_val_domain_boundary (tt(3)(dom), bound, res, pos_res, ncond(2)) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if ((*p_cmp[i])(0)==1)
						export_tau_val_domain_boundary (tt(1)(dom), bound, res, pos_res, ncond(i)) ;
					if ((*p_cmp[i])(0)==2)
						export_tau_val_domain_boundary (tt(2)(dom), bound, res, pos_res, ncond(i)) ;
					if ((*p_cmp[i])(0)==3)
						export_tau_val_domain_boundary (tt(3)(dom), bound, res, pos_res, ncond(i)) ;
				}
			break ;
		case 2 :
			if (tt.get_n_comp()==6) {
				if (n_cmp==-1) {
					export_tau_val_domain_boundary (tt(1,1)(dom), bound, res, pos_res, ncond(0)) ;
					export_tau_val_domain_boundary (tt(1,2)(dom), bound, res, pos_res, ncond(1)) ;
					export_tau_val_domain_boundary (tt(1,3)(dom), bound, res, pos_res, ncond(2)) ;
					export_tau_val_domain_boundary (tt(2,2)(dom), bound, res, pos_res, ncond(3)) ;
					export_tau_val_domain_boundary (tt(2,3)(dom), bound, res, pos_res, ncond(4)) ;
					export_tau_val_domain_boundary (tt(3,3)(dom), bound,  res, pos_res, ncond(5)) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain_boundary (tt(1, 1)(dom), bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain_boundary (tt(1, 2)(dom), bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(1, 3)(dom), bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain_boundary (tt(2, 2)(dom), bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(2, 3)(dom), bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(3, 3)(dom), bound, res, pos_res, ncond(i)) ;
				}
			}
			//  not symetric
			if (tt.get_n_comp()==9) {
				if (n_cmp==-1) {
					export_tau_val_domain_boundary (tt(1,1)(dom), bound, res, pos_res, ncond(0)) ;
					export_tau_val_domain_boundary (tt(1,2)(dom), bound, res, pos_res, ncond(1)) ;
					export_tau_val_domain_boundary (tt(1,3)(dom), bound, res, pos_res, ncond(2)) ;
					export_tau_val_domain_boundary (tt(2,1)(dom), bound, res, pos_res, ncond(3)) ;
					export_tau_val_domain_boundary (tt(2,2)(dom), bound, res, pos_res, ncond(4)) ;
					export_tau_val_domain_boundary (tt(2,3)(dom), bound,  res, pos_res, ncond(5)) ;
					export_tau_val_domain_boundary (tt(3,1)(dom), bound, res, pos_res, ncond(6)) ;
					export_tau_val_domain_boundary (tt(3,2)(dom), bound, res, pos_res, ncond(7)) ;
					export_tau_val_domain_boundary (tt(3,3)(dom), bound,  res, pos_res, ncond(8)) ;
					
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain_boundary (tt(1, 1)(dom), bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain_boundary (tt(1, 2)(dom), bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(1, 3)(dom), bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain_boundary (tt(2, 1)(dom), bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain_boundary (tt(2, 2)(dom), bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(2, 3)(dom), bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain_boundary (tt(3, 1)(dom), bound,  res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain_boundary (tt(3, 2)(dom), bound, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain_boundary (tt(3, 3)(dom), bound, res, pos_res, ncond(i)) ;
				}
			}
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_polar_periodic_nucleus::export_tau_boundary" << endl ;
			abort() ;
			break ;
	}
}}
