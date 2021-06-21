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
int Domain_polar_nucleus::nbr_conditions_val_domain (const Val_domain& so, int mquant, int llim, int order) const {
	int res = 0 ;
	
	Index pos (nbr_coefs) ;
	do {
		bool indic = true ;
		// Get base in theta :
		int baset = (*so.get_base().bases_1d[1]) (0) ;
		int lquant ;
		switch (baset) {
					case COS_EVEN:
						if ((pos(1)==0) && (mquant!=0))
							indic = false ;
						lquant = 2*pos(1) ;
						break ;
					case COS_ODD:
						if ((pos(1)==nbr_coefs(1)-1) || ((pos(1)==0) && (mquant!=0)))
							indic = false ;
						lquant = 2*pos(1)+1 ;
						break ;
					case SIN_EVEN:
						if (((pos(1)==1) && (mquant>1)) || (pos(1)==0) || (pos(1)==nbr_coefs(1)-1)) 
							indic = false  ;
						lquant = 2*pos(1) ;
						break ;
					case SIN_ODD:
						if  (((pos(1)==0) && (mquant>1)) || (pos(1)==nbr_coefs(1)-1))
							indic = false ;
						lquant = 2*pos(1)+1 ;
						break ;
					default:
						cerr << "Unknow theta basis in Domain_polar_nucleus::nbr_unknowns_val_domain" << endl ;
						abort() ;
		}

		int max = 0;
		if (indic) {
			// Base in r :
			int baser = (*so.get_base().bases_1d[0]) (pos(1)) ;
			 switch (baser) {
		case CHEB_EVEN : 
		    switch (baset) {
			case COS_EVEN :
			    if ((pos(0)==0) && ((lquant>llim) || (mquant!=0)))
				  indic = false ;
			    break ;
			case SIN_EVEN :
			    if (pos(0)==0)
				  indic = false ;
			    break ;
			default :
			  cerr << "Strange base in Domain_polar_nucleus::nbr_conditions_val_domain" << endl ;
			  abort() ;
		    }
		    max = nbr_coefs(0) ;
		    break ;
	    case LEG_EVEN : 
		    switch (baset) {
			case COS_EVEN :
			    if ((pos(0)==0) && ((lquant>llim) || (mquant!=0)))
				  indic = false ;
			    break ;
			case SIN_EVEN :
			    if (pos(0)==0)
				  indic = false ;
			    break ;
			default :
			  cerr << "Strange base in Domain_polar_nucleus::nbr_conditions_val_domain" << endl ;
			  abort() ;
		    }
		    max = nbr_coefs(0) ;
		    break ;
	    case CHEB_ODD :  
		    switch (baset) {
		      case SIN_ODD :
			if ((pos(0)==nbr_coefs(0)-1) || ((pos(0)==0) && (lquant>llim+1)))
			  indic = false ;
			break ;
		      case COS_ODD :
			if ((pos(0)==nbr_coefs(0)-1) || ((pos(0)==0) && (lquant>llim+1)))
			  indic = false ;
			break ;
		      default :
			  cerr << "Strange base in Domain_polar_nucleus::nbr_conditions_val_domain" << endl ;
			  abort() ;
		      }
		    max = nbr_coefs(0)-1 ;
		    break ;
	    case LEG_ODD :
		   switch (baset) {
		      case SIN_ODD :
			if ((pos(0)==nbr_coefs(0)-1) || ((pos(0)==0) && (lquant>llim+1)))
			  indic = false ;
			break ;
		      case COS_ODD :
			if ((pos(0)==nbr_coefs(0)-1) || ((pos(0)==0) && (lquant>llim+1)))
			  indic = false ;
			break ;
		      default :
			  cerr << "Strange base in Domain_polar_nucleus::nbr_conditions_val_domain" << endl ;
			  abort() ;
		      }
		    max = nbr_coefs(0)-1 ;
		    break ;
	    default :
	      cerr << "Unknown base in Domain_polar_nucleus::nbr_conditions_val_domain" << endl ;
	      abort() ;
	}
		}

		// Order with respect to r :
		int lim = 0 ;
		switch (order) {
		  case 2 :
		      lim = max-2 ;
		      break ;
		  case 1 :
		  {
		      if (pos(1)==0)
			lim = max-2 ;
		      else
			lim = max -1 ;}
		      break ;
		  case 0 :
		      lim = max-1 ;
		      break ;
		  default :
		      cerr << "Unknown case in Domain_polar_nucleus_nbr_conditions" << endl ;
		      abort() ;
		}
		
		if (pos(0)>lim)
			indic = false ;

		if (indic)
			res ++ ;
	}
	while (pos.inc()) ;

	return res ;
}

Array<int> Domain_polar_nucleus::nbr_conditions (const Tensor& tt, int dom, int order, int n_cmp, Array<int>**) const {

	int size = (n_cmp==-1) ? tt.get_n_comp() : n_cmp ;
	Array<int> res (size) ;
	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			if (!tt.is_m_quant_affected())  
			  res.set(0) = nbr_conditions_val_domain (tt()(dom), 0, 0, order) ;
			else
			  res.set(0) = nbr_conditions_val_domain (tt()(dom), tt.get_parameters()->get_m_quant(), 0, order) ;
			break ;
		default :
			cerr << "Valence " << val << " not implemented in Domain_polar_nucleus::nbr_conditions" << endl ;
			break ;
	}
	return res ;
}}
