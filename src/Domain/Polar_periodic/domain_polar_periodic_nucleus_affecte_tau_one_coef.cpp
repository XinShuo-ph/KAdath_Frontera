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
void Domain_polar_periodic_nucleus::affecte_tau_one_coef_val_domain (Val_domain& so, int llim, int cc, int& conte) const {

	int lquant ;

	so.is_zero = false ;
	so.allocate_coef() ;
	*so.cf=0. ;
	Index pos_cf(nbr_coefs) ;

	bool found = false ;

	// Positions of the Galerkin basis
	Index pos_gal_r (nbr_coefs) ;
	double  fact_r ;
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
			cerr << "Unknown basis in Domain_polar_periodic_nucleus_affecte_tau_val_domain" << endl ;
			abort() ;
	}
	// Loop on time
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
			cerr << "Unknown theta basis in Domain_polar_periodic_nucleus_affecte_tau_val_domain" << endl ;
			abort() ;
	}

	for (int j=minj ; j<maxj ; j++) {	
		int baser = (*so.get_base().bases_1d[0]) (j,k) ;
	
		// Compute lquant :
		switch (baset) {
		case COS_EVEN :
			lquant = 2*j ;
			break ;
		case COS_ODD :
			lquant = 2*j+1 ;
			break ;
		case SIN_EVEN :
			lquant = 2*j ;
			break ;
		case SIN_ODD :
			lquant = 2*j+1 ;
			break ;
		
		default :
			cerr << "Unknown theta basis in Domain_polar_periodic_nucleus_affecte_tau_val_domain" << endl ;
			abort() ;
	}

		pos_cf.set(1) = j ;

		int minr, maxr ;
		switch (baser) {
		case CHEB_EVEN :
			minr=(lquant<=llim) ? 0 : 1 ;
			maxr=nbr_coefs(0) ;
			break ;
		case CHEB_ODD :
			minr= (lquant<=llim) ? 0 : 1;
			maxr=nbr_coefs(0)-1 ;
			break ;
		case LEG_EVEN :
			minr=(lquant<=llim) ? 0 : 1 ;
			maxr=nbr_coefs(0) ;
			break ;
		case LEG_ODD :
			minr=(lquant<=llim) ? 0 : 1 ;
			maxr=nbr_coefs(0)-1 ;
			break ;
		default :
			cerr << "Unknown r basis in Domain_polar_periodic_nucleus_affecte_tau_val_domain" << endl ;
			abort() ;
		}

		// Loop on r :
		for (int i=minr ; i<maxr ; i++) {
			pos_cf.set(0) = i ;
			// No garlekin
			if (lquant<=llim) {
				if (conte==cc)  {
					found = true ;
					so.cf->set(pos_cf) = 1. ;
					}
				conte ++ ;
			}
			else {	if (conte==cc) {
				found = true ;
				// Galerkin en r
				pos_gal_r = pos_cf ;
				pos_gal_r.set(0) = 0 ;
				switch (baser) {
					case CHEB_EVEN :
						fact_r = - pow(-1, i) ;
						break ;
					case LEG_EVEN : {
						fact_r = -1. ;
						for (int t=0 ; t<i ; t++)
							fact_r *= -double(2*t+1)/double(2*t+2) ;
						}
						break ;
					 case CHEB_ODD :
						    fact_r = - (2*i+1) * pow(-1, i) ;
						    break ;
					case LEG_ODD : {
						    fact_r = -1. ;
						    for (int t=0 ; t<i ; t++)
						      fact_r *= -double(2*t+3)/double(2*t+2) ;
						    }
						    break ;

					default :
						cerr << "Strange base in Domain_polar_periodic_nucleus::affecte_tau_val_domain" << endl ;
						abort()  ;
				}
				so.cf->set(pos_cf) = 1 ;
				so.cf->set(pos_gal_r) += fact_r ;
				}
				conte ++ ;
			}
		 } // end loop i
		} // end loop j
		} // end loop k
	// If not found put to zero :
	if (!found)
		so.set_zero() ;
}


void Domain_polar_periodic_nucleus::affecte_tau_one_coef (Tensor& tt, int dom, int cc, int& pos_cf) const {

	// Check right domain
	assert (tt.get_space().get_domain(dom)==this) ;

	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			affecte_tau_one_coef_val_domain (tt.set().set_domain(dom), 0, cc, pos_cf) ;
			break ;
		case 1 :
			affecte_tau_one_coef_val_domain (tt.set(1).set_domain(dom), 2, cc, pos_cf) ;
			affecte_tau_one_coef_val_domain (tt.set(2).set_domain(dom), 2, cc, pos_cf) ;
			affecte_tau_one_coef_val_domain (tt.set(3).set_domain(dom), 2, cc, pos_cf) ;
			break ;
		case 2 :
			// symetric
			if (tt.get_n_comp()==6) {
				affecte_tau_one_coef_val_domain (tt.set(1,1).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,2).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,3).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,2).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,3).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,3).set_domain(dom), 2, cc, pos_cf) ;
			}
			//  not symetric
			if (tt.get_n_comp()==9) {
				affecte_tau_one_coef_val_domain (tt.set(1,1).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,2).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,3).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,1).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,2).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,3).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,1).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,2).set_domain(dom), 2, cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,3).set_domain(dom), 2, cc, pos_cf) ;
			}
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_polar_nucleus::affecte_tau" << endl ;
			abort() ;
			break ;
	}
}}
