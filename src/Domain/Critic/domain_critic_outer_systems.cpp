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
#include "point.hpp"
#include "array.hpp"
#include "val_domain.hpp"
namespace Kadath {
void Domain_critic_outer::find_other_dom (int dom, int bound, int& other_dom, int& other_bound) const {

	switch (bound) {
		case INNER_BC:
			other_dom = dom -1 ;
			other_bound = OUTER_BC ;
			break ;
		default:
			cerr << "Unknown boundary case in Domain_critic_outer::find_other_dom" << endl ;
			abort() ;
		}
}

double Domain_critic_outer::val_boundary (int bound, const Val_domain& so, const Index& pos_cf) const {

	if (so.check_if_zero())
		return 0. ;
	
	else {
	so.coef() ;
	double res = 0 ;
	Index copie_pos (pos_cf) ;

	switch (bound) {
		case INNER_BC :
			for (int i=0 ; i<nbr_coefs(0) ; i++) {
				copie_pos.set(0) = i ;
				if (i%2==0)
					res += so.get_coef(copie_pos) ;
				else
					res -= so.get_coef(copie_pos) ;
			}
			break ;
		case OUTER_BC : {
			for (int i=0 ; i<nbr_coefs(0) ; i++) {
				copie_pos.set(0) = i ;
				res += so.get_coef(copie_pos) ;
				}
			}
			break ;
		default :
			cerr << "Unknown boundary type in Domain_critic_outer::val_boundary" << endl ;
			abort() ;
	}
	return res ;
	}
}

int Domain_critic_outer::nbr_points_boundary (int bound, const Base_spectral& bb) const {

	if ((bound!=INNER_BC) && (bound!=OUTER_BC)) {
	  cerr << "Unknown boundary in Domain_critic_outer::nbr_points_boundary" << endl ;
	  abort() ;
	}

	int base_t = (*bb.bases_1d[1])(0) ;
	if ((base_t!=COSSIN_EVEN) && (base_t!=COSSIN_ODD)) {
		cerr << "Unknown base in Domain_critic_outer::nbr_conditions_val_domain_boundary" << endl ;
		abort() ;
	}
	return nbr_coefs(1)-2 ;
}

void Domain_critic_outer::do_which_points_boundary (int bound, const Base_spectral&, Index** which_coef, int start) const {

	int pos_which = start ;
	Index pos (nbr_points) ;

	switch (bound) {
		case INNER_BC :
			pos.set(0) = 0 ;
			break ;
		case OUTER_BC :
			pos.set(0) = nbr_points(0)-1 ;
			break ;
		default :
			cerr << "Unknown boundary in Domain_critic_outer::do_which_points_boundary" << endl ;
			abort() ;
	}

	for (int j=0 ; j<nbr_points(1) ; j++) {
		pos.set(1) = j ;
		which_coef[pos_which] = new Index(pos) ;
		pos_which ++ ;
	}
}}
