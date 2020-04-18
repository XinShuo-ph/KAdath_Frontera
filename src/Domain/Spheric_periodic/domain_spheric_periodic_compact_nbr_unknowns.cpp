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
#include "spheric_periodic.hpp"
#include "array_math.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "tensor.hpp"

namespace Kadath {
int Domain_spheric_periodic_compact::nbr_unknowns_val_domain (const Val_domain& so) const {
	
	int res = 0 ;
	int baset = (*so.get_base().bases_1d[1])(0) ;
	Index pos (nbr_coefs) ;
	do {
		bool indic = true ;
		switch (baset) {
			case COS_EVEN:
				break ;
			case COS_ODD:
				if (pos(1)==nbr_coefs(1)-1)
					indic = false ;
				break ;		
			case SIN_ODD:
				if (pos(1)==nbr_coefs(1)-1)
					indic = false ;
				break ;
			case SIN_EVEN:
				if ((pos(1)==0) || (pos(1)==nbr_coefs(1)-1))
					indic = false ;
				break ;
			case COS:
				break ;
			default:
				cerr << "Unknow time basis in Domain_spheric_periodic_compact::nbr_unknowns_val_domain" << endl ;
				abort() ;
		}
		if (indic)
		  res ++ ;
	}
	while (pos.inc()) ;

	return res ;
}

int Domain_spheric_periodic_compact::nbr_unknowns (const Tensor& tt, int dom) const {

	// Check right domain
	assert (tt.get_space().get_domain(dom)==this) ;

	int res = 0 ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			res += nbr_unknowns_val_domain (tt()(dom)) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_spheric_periodic_compact::nbr_unknowns" << endl ;
			break ;
	}
	return res ;
}}
