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

#include "system_of_eqs.hpp"
#include "ope_eq.hpp"
#include "term_eq.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
namespace Kadath {
//    const Domain* dom ;
//    int ndom ;
//    int n_ope ;
//    Memory_mapped_array<Ope_eq*> parts ;
//    bool called ;
//    int n_comp ;
//    int n_cond_tot ;
//    Array<int>* n_cond ;
//
//    int n_cmp_used ;
//    Memory_mapped_array<Array<int>*> p_cmp_used ;

Equation::Equation (const Domain* zedom, int nd, int np, int nused, Array<int>** pused) :
	dom{zedom}, ndom{nd}, n_ope{np}, parts{n_ope,initialize}, called{false}, n_cond{nullptr}, n_cond_tot{},
	n_cmp_used{nused}, p_cmp_used{pused}
{}


Equation::~Equation() {
    for(auto & p : parts) safe_delete(p);
	if (called) delete n_cond ;
}

void Equation::apply(int& conte, Term_eq** res) {
	int old_conte = conte ;
	for (int i=0 ; i<n_ope ; i++) {
		if (res[conte]!=nullptr)
			*res[conte] = parts[i]->action() ; //NOT THREAD SAFE
		else
			res[conte] = new Term_eq (parts[i]->action()) ; //NOT THREAD SAFE
		conte ++ ;
	}

	if (!called) {
		Tensor copie (res[old_conte]->get_val_t()) ;
		n_cond = new Array<int> (do_nbr_conditions(copie)) ;
		n_comp = n_cond->get_size(0) ;
		n_cond_tot = 0 ;
		for (int i=0 ; i<n_cond->get_size(0) ; i++)
			n_cond_tot += (*n_cond)(i) ;
		called = true ;
	}
}}
