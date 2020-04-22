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
#include "tensor_impl.hpp"
#include "tensor.hpp"

namespace Kadath {
void Domain_bispheric_chi_first::affecte_tau_one_coef_val_domain (Val_domain& so,  int cc, int& conte) const {

	int basep = (*so.get_base().bases_1d[2]) (0) ;

	so.is_zero = false ;
	so.allocate_coef() ;
	*so.cf=0. ;
	Index pos_cf(nbr_coefs) ;

	bool found = false ;
	do {	
		// Check if true :
		bool indic = true ;
		switch (basep) {
		case COS :
			// Last odd ones
			if ((pos_cf(2)%2==1) && (pos_cf(1)==nbr_coefs(1)-1))
				indic = false ;
			// Regularity for even ones :
			if ((pos_cf(2)!=0) && (pos_cf(2)%2==0) && (pos_cf(1)==0))
				indic = false ;
			break ;
		case SIN :
			// sin(0)
			if ((pos_cf(2)==0) || (pos_cf(2)==nbr_coefs(2)-1))
				indic = false ;
			// Last odd ones :
			if ((pos_cf(2)%2==1) && (pos_cf(1)==nbr_coefs(1)-1))
				indic = false ;
			// Regularity for even ones :
			if ((pos_cf(2)%2==0) && (pos_cf(1)==0))
				indic = false ;
			break ;
		default :
			cerr << "Unknown phi basis in Domain_bispheric_chi_first::affecte_tau_val_domain" << endl ;
			abort() ;
		}

		if (indic) {
			if (conte==cc) {
				found = true ;
				so.cf->set(pos_cf) = 1;
				// Regularity on the axis :
				if ((pos_cf(2)%2==0) && (pos_cf(2)!=0)) {
					Index pos_galerkin (pos_cf) ;
					pos_galerkin.set(1) = 0 ;
					double valreg ;
					int basechi = (*so.get_base().bases_1d[1])(pos_cf(2)) ;
					switch (basechi) {
						case CHEB_EVEN :
							valreg = (pos_cf(1)%2==0) ? -1 : 1  ;
							break ;
						case LEG_EVEN :
							valreg = 0.5 ;
							for (int i=1 ; i<pos_cf(1) ; i++)
								valreg *= - double(2*i+1)/double(2*i+2) ;
							break ;
						default :
							cerr << "Unknown base in Domain_bispheric_chi_first::affecte_one_coef" << endl ;
							abort() ;
						}
					so.cf->set(pos_galerkin) = valreg ;
					}
				}
			else
				so.cf->set(pos_cf) = 0. ;
			conte ++ ;
		}
	}
	while (pos_cf.inc()) ;

	// If not found put to zero :
	if (!found)
		so.set_zero() ;
}

void Domain_bispheric_chi_first::affecte_tau_one_coef (Tensor& tt, int dom, int cc, int& pos_cf) const {

	// Check right domain
	assert (tt.get_space().get_domain(dom)==this) ;

	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			affecte_tau_one_coef_val_domain (tt.set().set_domain(dom), cc, pos_cf) ;
			break ;
		case 1 : {
			bool found = false ;
			// Cartesian basis
			if (tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) {
				affecte_tau_one_coef_val_domain (tt.set(1).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3).set_domain(dom), cc, pos_cf) ;
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of vector Domain_bispheric_chi_first::affecte_tau_one_coef" << endl ;
				abort() ;
			}
		}
			break ;
		case 2 : {
			bool found = false ;
			// Cartesian basis and symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==6)) {
				affecte_tau_one_coef_val_domain (tt.set(1,1).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,3).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,3).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,3).set_domain(dom), cc, pos_cf) ;
				found = true ;
			}
			// Cartesian basis and not symetric
			if ((tt.get_basis().get_basis(dom)==CARTESIAN_BASIS) && (tt.get_n_comp()==9)) {
				affecte_tau_one_coef_val_domain (tt.set(1,1).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(1,3).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,1).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(2,3).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,1).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,2).set_domain(dom), cc, pos_cf) ;
				affecte_tau_one_coef_val_domain (tt.set(3,3).set_domain(dom), cc, pos_cf) ;
				found = true ;
			}
			if (!found) {
				cerr << "Unknown type of 2-tensor Domain_bispheric_chi_first::affecte_tau_one_coef" << endl ;
				abort() ;
			}
		}
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_bispheric_chi_first::affecte_tau" << endl ;
			break ;
	}
}
}
