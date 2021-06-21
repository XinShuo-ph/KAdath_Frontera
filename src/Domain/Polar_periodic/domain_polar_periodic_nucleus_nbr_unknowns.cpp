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
#include "polar_periodic.hpp"
#include "point.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"

namespace Kadath {
int Domain_polar_periodic_nucleus::nbr_unknowns_val_domain (const Val_domain& so, int llim) const {
	
	int res = 0 ;
	Index pos (nbr_coefs) ;
	do {
		bool indic = true ;

		// Base in time
		int basetime = (*so.get_base().bases_1d[2]) (0) ;
		switch (basetime) {
					case COS:
						break ;
					case SIN:
						if ((pos(2)==0) || (pos(2)==nbr_coefs(2)-1))
							indic = false ;
						break ;
					default:
						cerr << "Unknow time basis in Domain_polar_polar_nucleus::nbr_unknowns_val_domain" << endl ;
						abort() ;
		}

		// Get base in theta :
		int baset = (*so.get_base().bases_1d[1]) (pos(2)) ;
		int lquant ;
		switch (baset) {
					case COS_EVEN:
						lquant = 2*pos(1) ;
						break ;
					case COS_ODD:
						if (pos(1)==nbr_coefs(1)-1)
							indic = false ;
						lquant = 2*pos(1)+1 ;
						break ;
					case SIN_EVEN:
						if ((pos(1)==0) || (pos(1)==nbr_coefs(1)-1)) 
							indic = false  ;
						lquant = 2*pos(1) ;
						break ;
					case SIN_ODD:
						if (pos(1)==nbr_coefs(1)-1)
							indic = false ;
						lquant = 2*pos(1)+1 ;
						break ;
					default:
						cerr << "Unknow theta basis in Domain_polar_periodic_nucleus::nbr_unknowns_val_domain" << endl ;
						abort() ;
		}
	
		if (indic) {		

		// Base in r :
		int baser = (*so.get_base().bases_1d[0]) (pos(1), pos(2)) ;
	        switch (baser) {
		case CHEB_EVEN : 
		    if ((pos(0)==0) && (lquant>llim))
				  indic = false ;
		    break ;
	    	case LEG_EVEN : 
		   if ((pos(0)==0) && (lquant>llim))
				  indic = false ;
		    break ;
	    case CHEB_ODD :  
		     if (((pos(0)==0) && (lquant>llim)) || (pos(0) == nbr_coefs(0)-1))
				  indic = false ;
		    break ;
	    case LEG_ODD :
		   if (((pos(0)==0) && (lquant>llim)) || (pos(0) == nbr_coefs(0)-1))
				  indic = false ;
		    break ;
	    default :
	      cerr << "Unknown base in Domain_polar_periodic_nucleus::nbr_unknowns_val_domain" << endl ;
	      abort() ;
	}
      }

		if (indic)
			res ++ ;
	}
	while (pos.inc()) ;

	return res ;
}

int Domain_polar_periodic_nucleus::nbr_unknowns (const Tensor& tt, int dom) const {

	// Check right domain
	assert (tt.get_space().get_domain(dom)==this) ;

	int res = 0 ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			res += nbr_unknowns_val_domain (tt()(dom), 0) ;
			break ;
		case 1 :
			res += nbr_unknowns_val_domain (tt(1)(dom), 2) ;
			res += nbr_unknowns_val_domain (tt(2)(dom), 2) ;
			res += nbr_unknowns_val_domain (tt(3)(dom), 2) ;
			break ;
		case 2 :
			if (tt.get_n_comp()==6) {
				res += nbr_unknowns_val_domain (tt(1,1)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(1,2)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(1,3)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(2,2)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(2,3)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(3,3)(dom), 2) ;
			}
			else {
				res += nbr_unknowns_val_domain (tt(1,1)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(1,2)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(1,3)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(2,1)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(2,2)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(2,3)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(3,1)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(3,2)(dom), 2) ;
				res += nbr_unknowns_val_domain (tt(3,3)(dom), 2) ;
			}
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_polar_periodic_nucleus::nbr_unknowns" << endl ;
			break ;
	}
	return res ;
}}
