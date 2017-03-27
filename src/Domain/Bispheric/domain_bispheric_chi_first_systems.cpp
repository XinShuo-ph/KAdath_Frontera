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
#include "param.hpp"
#include "val_domain.hpp"

namespace Kadath {
void Domain_bispheric_chi_first::find_other_dom (int dom, int bound, int& other_dom, int& other_bound) const {
  
	switch (bound) {
		case CHI_ONE_BC :
			if (eta_lim<0)
				other_dom = dom+1 ;
			else
				other_dom = dom-1 ;
			other_bound = CHI_ONE_BC ;
			break ;
		default:
			cerr << "Unknown boundary case in Domain_bispheric_chi_first::find_other_dom" << endl ;
			abort() ;
		}
}

double Domain_bispheric_chi_first::val_boundary (int bound, const Val_domain& so, const Index& pos_cf) const {

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
		case OUTER_BC :
			for (int i=0 ; i<nbr_coefs(0) ; i++) {
				copie_pos.set(0) = i ;
				res += so.get_coef(copie_pos) ;
				}
			break ;
		case CHI_ONE_BC :
			for (int j=0 ; j<nbr_coefs(1) ; j++) {
				copie_pos.set(1) = j ;
				res += so.get_coef(copie_pos) ;
			}
			break ;
		default :
			cerr << "Unknown boundary type in Domain_bispheric_chi_first::val_boundary" << endl ;
			abort() ;
	}
	return res ;
	}
}

int Domain_bispheric_chi_first::nbr_points_boundary (int bound, const Base_spectral& bb) const {
	int basep = (*bb.bases_1d[2]) (0) ;
	int nbrphi = (basep==COS) ? nbr_points(2) : nbr_points(2)-2 ;
	int res ;
	switch (bound) {
		case INNER_BC :
			res = nbrphi * (nbr_points(1)-1)  ;
			break ;
		case OUTER_BC : 
			res = nbrphi * (nbr_points(1) -1) ;
			break ;
		default :
			cerr << "Domain_bispheric_chi_first::nbr_points_boundary not yet implemented for boundary " << bound << endl ;
			abort() ;
		}
	if (basep==COS)
		res ++ ;

	return res ;
}

void Domain_bispheric_chi_first::do_which_points_boundary (int bound, const Base_spectral& bb, Index** which_coef, int start) const {
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
			cerr << "Unknown boundary in Domain_bispheric_chi_first::do_which_points_inside" << endl ;
			abort() ;
	}

	//Look at the symetrie :
	int basep = (*bb.bases_1d[2]) (0) ;
	int mink = (basep==COS) ? 0 : 1 ;
	int maxk = (basep==COS) ? nbr_points(2) : nbr_points(2)-1 ;

	// Point on the axis :
	if (basep==COS) {
		pos.set(2) = 0 ;
		pos.set(1) = 0 ;
		which_coef[pos_which] = new Index(pos) ;
		pos_which ++ ;
	}

	for (int k=mink ; k<maxk ; k++) {
		pos.set(2) = k ;
		for (int j=1 ; j<nbr_points(1) ; j++) {
			pos.set(1) = j ;
			which_coef[pos_which] = new Index(pos) ;
			pos_which ++ ;
		}
	}
}
}
