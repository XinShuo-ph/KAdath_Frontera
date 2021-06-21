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
#include "point.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
void Domain_polar_shell::export_tau_val_domain (const Val_domain& so, int mquant, int order, Array<double>& sec, int& pos_sec, int ncond) const {

	if (so.check_if_zero()) 
		pos_sec += ncond ;
	else {

	so.coef() ;
	Index pos_cf (nbr_coefs) ;
	Index pos_galerkin (nbr_coefs) ;

	// Loop on theta
	int baset = (*so.get_base().bases_1d[1]) (0) ;
	for (int j=0 ; j<nbr_coefs(1) ; j++) {
		pos_cf.set(1) = j ;
		// Loop on r :
		for (int i=0 ; i<nbr_coefs(0)-order ; i++) {
			pos_cf.set(0) = i ;
			switch (baset) {
				case COS_EVEN:
					if (mquant==0) {
						sec.set(pos_sec) = (*so.cf)(pos_cf) ;
						pos_sec ++ ;
						}
					else if (j!=0) {
						// Galerkin base
						pos_galerkin = pos_cf ;
						pos_galerkin.set(1) = 0 ;
						sec.set(pos_sec) = (*so.cf)(pos_cf) 
								-2*(*so.cf)(pos_galerkin) ;
						pos_sec ++ ;
						}	
					break ;
				case COS_ODD:
					if (j!=nbr_coefs(1)-1) {
						if (mquant==0) {
							sec.set(pos_sec) = (*so.cf)(pos_cf) ;
							pos_sec ++ ;
						}
						else if (j!=0) {
							// Galerkin base
							pos_galerkin = pos_cf ;
							pos_galerkin.set(1) = 0 ;
							sec.set(pos_sec) = (*so.cf)(pos_cf) 
								-(*so.cf)(pos_galerkin) ;
							pos_sec ++ ;
						}
					}
					break ;
				case SIN_EVEN:
					if ((j!=0) && (j!=nbr_coefs(1)-1)) {
						if (mquant<=1) {
							sec.set(pos_sec) = (*so.cf)(pos_cf) ;
							pos_sec ++ ;
						}
						else if (j!=1) {
							// Galerkin base
							pos_galerkin = pos_cf ;
							pos_galerkin.set(1) = 1 ;
							sec.set(pos_sec) = (*so.cf)(pos_cf) 
								- j*(*so.cf)(pos_galerkin) ;
							pos_sec ++ ;
						}
					}
					break ;
				case SIN_ODD:
					if (j!=nbr_coefs(1)-1) {
						if (mquant<=1) {
						sec.set(pos_sec) = (*so.cf)(pos_cf) ;
							pos_sec ++ ;
						}
						else if (j!=0) {
							// Galerkin base
							pos_galerkin = pos_cf ;
							pos_galerkin.set(1) = 0 ;
							sec.set(pos_sec) = (*so.cf)(pos_cf) 
								- (2*j+1)*(*so.cf)(pos_galerkin) ;
							pos_sec ++ ;
						}
					}
					break ;
				default:
					cerr << "Unknow theta basis in Domain_polar_shell::export_tau_val_domain" << endl ;
					abort() ;
				}
			}
		}
	}
}

void Domain_polar_shell::export_tau (const Tensor& tt, int dom, int order, Array<double>& res, int& pos_res, const Array<int>& ncond,
										int, Array<int>**) const {
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			if (!tt.is_m_quant_affected())
			  export_tau_val_domain (tt()(dom), 0, order, res, pos_res, ncond(0)) ;
			else 
			    export_tau_val_domain (tt()(dom), tt.get_parameters()->get_m_quant(), order, res, pos_res, ncond(0)) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_polar_shell::export_tau" << endl ;
			break ;
	}
}}
