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
int Domain_polar_compact::nbr_conditions_val_domain (const Val_domain& so, int mquant, int order) const {
	int res = 0 ;
	int rlim = order ;
	Index pos (nbr_coefs) ;
	do {
	  
	
		bool indic = true ;
		// Get base in theta :
		int baset = (*so.get_base().bases_1d[1]) (0) ;	
	//	if ((pos(1)!=0) && (order==1) && (baset==COS_EVEN)) 
		//    rlim = 2 ;
		
		switch (baset) {
					case COS_EVEN:
						if ((pos(1)==0) && (mquant!=0))
							indic = false ;
						break ;
					case COS_ODD:
						if ((pos(1)==nbr_coefs(1)-1) || ((pos(1)==0) && (mquant!=0)))
							indic = false ;
						break ;
					case SIN_EVEN:
						 if (((pos(1)==1) && (mquant>1)) || (pos(1)==0) || (pos(1)==nbr_coefs(1)-1))
							indic = false  ;
						break ;
					case SIN_ODD:
						if (((pos(1)==0) && (mquant>1)) || (pos(1)==nbr_coefs(1)-1))
							indic = false ;
						break ;
					default:
						cerr << "Unknow theta basis in Domain_polar_compact::nbr_conditions_val_domain" << endl ;
						abort() ;
		}
		// Order with respect to r :
		if (pos(0)>nbr_coefs(0)-rlim-1)
			indic = false ;

		if (indic)
			res ++ ;
	}
	while (pos.inc()) ;

	return res ;
}

Array<int> Domain_polar_compact::nbr_conditions (const Tensor& tt, int dom, int order, int n_cmp, Array<int>**) const {

	int size = (n_cmp==-1) ? tt.get_n_comp() : n_cmp ;
	Array<int> res (size) ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			if (!tt.is_m_quant_affected())
			  res.set(0) = nbr_conditions_val_domain (tt()(dom), 0, order) ;
			else 
			  res.set(0) = nbr_conditions_val_domain (tt()(dom), tt.get_parameters()->get_m_quant(), order) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_polar_compact::nbr_conditions" << endl ;
			break ;
	}
	return res ;
}}
