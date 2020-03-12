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
#include "utilities.hpp"
#include "homothetic.hpp"
#include "array_math.hpp"
#include "val_domain.hpp"
#include "scalar.hpp"
#include "vector.hpp"
namespace Kadath {
Term_eq Domain_shell_inner_homothetic::integ_term_eq (const Term_eq& so, int bound) const {
    // Check it is a tensor
	if (so.get_type_data() != TERM_T) {
		cerr << "integ_term_eq only defined with respect for a tensor" << endl ;
		abort() ;
	}

	if (so.get_val_t().get_n_comp() != 1) {
		cerr << "integ_term_eq only defined with respect to a scalar" << endl ;
		abort() ;
	}
	
	assert (so.get_dom()==num_dom) ;
	

	Term_eq rrso (mult_r_term_eq(mult_r_term_eq(so))) ;
	
	
	
	double resval = 0 ;
	{// The value
	Array<int> ind (so.get_val_t().indices(0)) ;
	Val_domain value (mult_sin_theta(rrso.get_val_t()(ind)(num_dom))) ;
	if (value.check_if_zero()) 
		resval= 0. ;
	else  {

	  int baset = (*value.get_base().bases_1d[1]) (0) ;
	  Index pcf (nbr_coefs) ;
	switch (baset) {
	  case COS_ODD :
	break ;
	case SIN_EVEN :
	  break ;
	case COS_EVEN : {
	resval += M_PI*val_boundary(bound, value, pcf) ;
	break ;
    }
    case SIN_ODD : {
	  for (int j=0 ; j<nbr_coefs(1) ; j++) {
	    pcf.set(1) = j ;
	    resval += 2./(2*double(j)+1) * val_boundary(bound, value, pcf) ;
	  }
      break ;
    }
    
    default : 
      cerr << "Case not yet implemented in Domain_shell_inner_homothetic::integrale" << endl ;
      abort() ;
  }
 
	}
    resval *= 2*M_PI ;
	}
	
	if (rrso.get_p_der_t()!=0x0) {  
	// The value
	Array<int> ind (so.get_val_t().indices(0)) ;
	Val_domain valueder (mult_sin_theta(rrso.get_der_t()(ind)(num_dom))) ;
	double resder = 0;
	if (valueder.check_if_zero()) 
		resder= 0. ;
	else  {

	  int baset = (*valueder.get_base().bases_1d[1]) (0) ;
	  Index pcf (nbr_coefs) ;
	switch (baset) {
	  case COS_ODD :
	break ;
	case SIN_EVEN :
	  break ;
	case COS_EVEN : {
	resder += M_PI*val_boundary(bound, valueder, pcf) ;
	break ;
    }
    case SIN_ODD : {
	  for (int j=0 ; j<nbr_coefs(1) ; j++) {
	    pcf.set(1) = j ;
	    resder += 2./(2*double(j)+1) * val_boundary(bound, valueder, pcf) ;
	  }
      break ;
    }
    default : 
      cerr << "Case not yet implemented in Domain_shell_inner_homothetic::integrale" << endl ;
      abort() ;
  }
	resder *= 2*M_PI ;
	}
	  return Term_eq (num_dom, resval, resder) ;
	}
	else {
	  return Term_eq (num_dom, resval) ;
	}
}

}
