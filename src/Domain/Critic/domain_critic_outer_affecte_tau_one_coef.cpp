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
#include "critic.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"
namespace Kadath {
void Domain_critic_outer::affecte_tau_one_coef_val_domain (Val_domain& so, int cc, int& conte) const {

	so.is_zero = false ;
	so.allocate_coef() ;
	*so.cf=0. ;
	Index pos(nbr_coefs) ;

	bool found = false ;

	do  {
		bool indic = true ;
		int base_t = (*so.get_base().bases_1d[1])(0) ;

		assert ((base_t==COSSIN_EVEN) || (base_t==COSSIN_ODD)) ;
		if (pos(1)==nbr_coefs(1)-1)
			indic = false ;
		if ((pos(1)==1) && (base_t==COSSIN_EVEN))
			indic = false ;
		if ((pos(1)==nbr_coefs(1)-2) && (base_t==COSSIN_ODD))
			indic = false ;

		if (indic) {
			if (conte==cc) {
				so.cf->set(pos) = 1;
				found = true ;
				}
			else
				so.cf->set(pos) = 0. ;
			conte ++ ;
		}
	}
	while (pos.inc()) ;

	// If not found put to zero :
	if (!found)
		so.set_zero() ;
}

void Domain_critic_outer::affecte_tau_one_coef (Tensor& tt, int dom, int cc, int& pos_cf) const {

	// Check right domain
	assert (tt.get_space().get_domain(dom)==this) ;

	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			affecte_tau_one_coef_val_domain (tt.set().set_domain(dom), cc, pos_cf) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_critic_outer::affecte_tau" << endl ;
			break ;
	}
}}
