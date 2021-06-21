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
int Domain_critic_outer::nbr_conditions_val_domain (const Val_domain& so, int order) const {

	int res = 0 ;

	int base_t = (*so.get_base().bases_1d[1])(0) ;
	assert ((base_t==COSSIN_EVEN) || (base_t==COSSIN_ODD)) ;

	//Loop on t
	int max = (base_t==COSSIN_EVEN) ? nbr_coefs(1)-1 : nbr_coefs(1)-2 ;
	for (int j=0 ; j<max ; j++)
		if ((j!=1) || (base_t!=COSSIN_EVEN))
		  res += nbr_coefs(0)-order ;
	return res ;
}

Array<int> Domain_critic_outer::nbr_conditions (const Tensor& tt, int dom, int order, int n_cmp, Array<int>**) const {

	int size = (n_cmp==-1) ? tt.get_n_comp() : n_cmp ;
	Array<int> res (size) ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			res.set(0) = nbr_conditions_val_domain (tt()(dom), order) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_shell::nbr_conditions" << endl ;
			break ;
	}
	return res ;
}}
