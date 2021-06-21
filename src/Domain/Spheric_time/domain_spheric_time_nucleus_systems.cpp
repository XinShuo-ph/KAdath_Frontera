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
#include "spheric_time.hpp"
#include "array.hpp"
#include "val_domain.hpp"
namespace Kadath {
void Domain_spheric_time_nucleus::find_other_dom (int dom, int bound, int& other_dom, int& other_bound) const {

	switch (bound) {
		case OUTER_BC:
			other_dom = dom +1 ;
			other_bound = INNER_BC ;
			break ;
		default:
			cerr << "Unknown boundary case in Domain_spheric_time_nucleus::find_other_dom" << endl ;
			abort() ;
		}
}

double Domain_spheric_time_nucleus::val_boundary (int bound, const Val_domain& so, const Index& pos_cf) const {

	if (so.check_if_zero())
		return 0. ;
	
	else {
	so.coef() ;
	double res = 0 ;
	Index copie_pos (pos_cf) ;
	switch (bound) {
		case OUTER_BC :
			for (int i=0 ; i<nbr_coefs(0) ; i++) {
				copie_pos.set(0) = i ;
				res += so.get_coef(copie_pos) ;
				}
			break ;
		case INNER_BC :
			for (int i=0 ; i<nbr_coefs(0) ; i++) {
				copie_pos.set(0) = i ;
				if (i%2==0)
				  res += so.get_coef(copie_pos) ;
				else
				  res -= so.get_coef(copie_pos) ;
				}
			break ;
		case TIME_INIT :
		  for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    copie_pos.set(1) = j ;
		    if (j%2==0)
			res += so.get_coef(copie_pos) ;
		    else
			res -= so.get_coef(copie_pos) ;
		}
		break ;
		default :
			cerr << "Unknown boundary type in Domain_spheric_time_nucleus::val_boundary" << endl ;
			abort() ;
	}

	return res ;
	}
}}
