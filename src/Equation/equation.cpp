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
namespace Kadath {
Equation::Equation (const Domain* zedom, int nd, int np, int nused, Array<int>** pused) : 
	dom(zedom), ndom(nd), 
	n_ope(np), called(false), n_cond(0x0), n_cmp_used(nused), p_cmp_used(pused) {

	parts = new Ope_eq*[n_ope] ;
	for (int i=0 ; i<n_ope ; i++)
		parts[i] = 0x0 ;
}

Equation::~Equation() {
	for (int i=0 ; i<n_ope ; i++)
		if (parts[i] !=0x0)
			delete parts[i] ;
	delete [] parts ;

	if (called) 	
		delete n_cond ;
}

void Equation::apply(int& conte, Term_eq** res) {
	int old_conte = conte ;
	for (int i=0 ; i<n_ope ; i++) {
		if (res[conte]!=0x0)
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
