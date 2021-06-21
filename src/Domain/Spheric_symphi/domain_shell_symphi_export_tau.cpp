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
#include "spheric_symphi.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
void Domain_shell_symphi::export_tau_val_domain (const Val_domain& so, int order, Array<double>& sec, int& pos_sec, int ncond) const {

	if (so.check_if_zero()) 
		pos_sec += ncond ;
	else {
	  
	so.coef() ;
	
	Index pos_cf (nbr_coefs) ;
	// Positions of the Galerkin basis
	Index pos_gal_t (nbr_coefs) ;
	double fact_t ;
	

	int mquant, kmin, kmax ;
		// Base in phi 
		int basep = (*so.get_base().bases_1d[2]) (0) ;
		switch (basep) {
					case COS_EVEN:
						kmin = 0 ;
						kmax = nbr_coefs(2)-1 ;
						break ;
					case COS_ODD:
						kmin = 0 ;
						kmax = nbr_coefs(2)-2 ;
						break ;
					case SIN_EVEN:
						kmin = 1 ;
						kmax = nbr_coefs(2)-2 ;
						break ;
					case SIN_ODD:
						kmin = 0 ;
						kmax = nbr_coefs(2)-2 ;
						break ;
					default:
						cerr << "Unknow phi basis in Domain_shell_symphi::export_tau_val_domain" << endl ;
						abort() ;
		}



	// Loop on phi :
	for (int k=kmin ; k<=kmax ; k++) {

			switch (basep) {
					case COS_EVEN:
						mquant = 2*k ;
						break ;
					case COS_ODD:
						mquant = 2*k+1 ;
						break ;
					case SIN_EVEN:
						mquant = 2*k ;
						break ;
					case SIN_ODD:
						mquant = 2*k+1 ;
						break ;
					default:
						cerr << "Unknow phi basis in Domain_shell_symphi::export_tau_val_domain" << endl ;
						abort() ;
			}



			pos_cf.set(2) = k ;
			// Loop on theta
			int baset = (*so.get_base().bases_1d[1]) (k) ;
			for (int j=0 ; j<nbr_coefs(1) ; j++) {	
				int baser = (*so.get_base().bases_1d[0]) (j, k) ;
				pos_cf.set(1) = j ;
				// Loop on r :
				for (int i=0 ; i<nbr_coefs(0)-order ; i++) {
					pos_cf.set(0) = i ;
					switch (baset) {
						case COS_EVEN :
							// No galerkin 
							if (mquant==0) {
								sec.set(pos_sec) = (*so.cf)(pos_cf) ;
								pos_sec ++ ;
							}
							else if (j!=0) {
							    // Need to use  Galerkin basis 
							    pos_gal_t = pos_cf ;
							    pos_gal_t.set(1) = 0 ;
							    fact_t = -2. ;
							    sec.set(pos_sec) = (*so.cf)(pos_cf) + fact_t*(*so.cf)(pos_gal_t) ;
							    pos_sec++ ;
							    }
						break ;
					case COS_ODD:
						// True coefs ?
						if (j!=nbr_coefs(1)-1)  {
							if (mquant==0) {
							  // No Galerkin :
							  sec.set(pos_sec) = (*so.cf)(pos_cf) ;
							  pos_sec ++ ;
							}
						      else if (j!=0) {
							    // Need to use  Galerkin basis 
							    pos_gal_t = pos_cf ;
							    pos_gal_t.set(1) = 0 ;
							    fact_t = -1. ;
							    sec.set(pos_sec) = (*so.cf)(pos_cf) + fact_t*(*so.cf)(pos_gal_t)  ;
							    pos_sec++ ;
							}
						}
						break ;
					case SIN_EVEN:
						if ((j!=0) && (j!=nbr_coefs(1)-1)) {
							if (mquant<=1) {
							      // No Galerkin
							      sec.set(pos_sec) = (*so.cf)(pos_cf) ;
							      pos_sec ++ ;
							} 
							else if (j!=1) {
								//  Galerkin
							    pos_gal_t = pos_cf ;
							    pos_gal_t.set(1) = 1 ;
							    fact_t = -j ;
							     sec.set(pos_sec) = (*so.cf)(pos_cf) + fact_t*(*so.cf)(pos_gal_t) ;
							    pos_sec++ ;
						 }
						}
						break ;
					case SIN_ODD:
						// True coefs ?
						if (j!=nbr_coefs(1)-1)  {
							if (mquant<=1)  {
							  // No Galerkin :
							  sec.set(pos_sec) = (*so.cf)(pos_cf) ;
							  pos_sec ++ ;
							}
						      else if (j!=0) {
							    // Need to use Galerkin basis
							    pos_gal_t = pos_cf ;
							    pos_gal_t.set(1) = 0 ;
							    fact_t = -(2*j+1) ;
							    sec.set(pos_sec) = (*so.cf)(pos_cf) + fact_t*(*so.cf)(pos_gal_t) ;
							    pos_sec++ ;
							}
						}
						break ;
					default:
						cerr << "Unknow theta basis in Domain_shell_symphi::export_tau_val_domain" << endl ;
						abort() ;
						}
					}
			}
		}
	}
}

void Domain_shell_symphi::export_tau (const Tensor& tt, int dom, int order, Array<double>& res, int& pos_res, const Array<int>& ncond,
										int n_cmp, Array<int>** p_cmp) const {
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
		      export_tau_val_domain (tt()(dom), order, res, pos_res, ncond(0)) ;
			break ;
		case 1 : {
			bool found = false ;
			// Cartesian basis
			if (tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) {
				if (n_cmp==-1) {
					export_tau_val_domain (tt(1)(dom), order, res, pos_res, ncond(0)) ;
					export_tau_val_domain (tt(2)(dom), order, res, pos_res, ncond(1)) ;
					export_tau_val_domain (tt(3)(dom), order, res, pos_res, ncond(2)) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if ((*p_cmp[i])(0)==1)
						export_tau_val_domain (tt(1)(dom), order, res, pos_res, ncond(i)) ;
					if ((*p_cmp[i])(0)==2)
						export_tau_val_domain (tt(2)(dom), order, res, pos_res, ncond(i)) ;
					if ((*p_cmp[i])(0)==3)
						export_tau_val_domain (tt(3)(dom), order, res, pos_res, ncond(i)) ;
				}
				found = true ;
			}
			
			if (!found) {
				cerr << "Unknown type of vector Domain_shell_symphi::export_tau" << endl ;
				abort() ;
			}
		}
			break ;
		case 2 : {
			bool found = false ;
			// Cartesian basis and symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==6)) {
				if (n_cmp==-1) {
					export_tau_val_domain (tt(1,1)(dom), order, res, pos_res, ncond(0)) ;
					export_tau_val_domain (tt(1,2)(dom), order, res, pos_res, ncond(1)) ;
					export_tau_val_domain (tt(1,3)(dom), order, res, pos_res, ncond(2)) ;
					export_tau_val_domain (tt(2,2)(dom), order, res, pos_res, ncond(3)) ;
					export_tau_val_domain (tt(2,3)(dom), order, res, pos_res, ncond(4)) ;
					export_tau_val_domain (tt(3,3)(dom), order, res, pos_res, ncond(5)) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain (tt(1, 1)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain (tt(1, 2)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain (tt(1, 3)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain (tt(2, 2)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain (tt(2, 3)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain (tt(3, 3)(dom), order, res, pos_res, ncond(i)) ;
				}
				found = true ;
			}
			// Cartesian basis and not symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==9)) {
				if (n_cmp==-1) {
					export_tau_val_domain (tt(1,1)(dom), order, res, pos_res, ncond(0)) ;
					export_tau_val_domain (tt(1,2)(dom), order, res, pos_res, ncond(1)) ;
					export_tau_val_domain (tt(1,3)(dom), order, res, pos_res, ncond(2)) ;
					export_tau_val_domain (tt(2,1)(dom), order, res, pos_res, ncond(3)) ;
					export_tau_val_domain (tt(2,2)(dom), order, res, pos_res, ncond(4)) ;
					export_tau_val_domain (tt(2,3)(dom), order, res, pos_res, ncond(5)) ;
					export_tau_val_domain (tt(3,1)(dom), order, res, pos_res, ncond(6)) ;
					export_tau_val_domain (tt(3,2)(dom), order, res, pos_res, ncond(7)) ;
					export_tau_val_domain (tt(3,3)(dom), order, res, pos_res, ncond(8)) ;
				}
				else for (int i=0 ; i<n_cmp ; i++) {
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain (tt(1, 1)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain (tt(1, 2)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==1) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain (tt(1, 3)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain (tt(2, 1)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain (tt(2, 2)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==2) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain (tt(2, 3)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==1))
						export_tau_val_domain (tt(3, 1)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==2))
						export_tau_val_domain (tt(3, 2)(dom), order, res, pos_res, ncond(i)) ;
					if (((*p_cmp[i])(0)==3) && ((*p_cmp[i])(1)==3))
						export_tau_val_domain (tt(3, 3)(dom), order, res, pos_res, ncond(i)) ;
				}
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of 2-tensor Domain_shell_symphi::export_tau" << endl ;
				abort() ;
			}
		}
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_shell_symphi::export_tau" << endl ;
			break ;
	}
}
}
