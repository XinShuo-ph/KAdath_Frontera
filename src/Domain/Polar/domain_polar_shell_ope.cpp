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
#include "array.hpp"
#include "polar.hpp"
#include "val_domain.hpp"
#include "array_math.hpp"
namespace Kadath {
int mult_cos_1d (int, Array<double>&) ;
int mult_sin_1d (int, Array<double>&) ;
int div_sin_1d (int, Array<double>&) ;
int mult_x_1d (int, Array<double>&) ;
int div_xp1_1d (int, Array<double>&) ;

Val_domain Domain_polar_shell::mult_cos_theta (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base = so.base ;	
	
	res.cf = (so.base.ope_1d(mult_cos_1d, 1, so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}

Val_domain Domain_polar_shell::mult_sin_theta (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;
	
	res.cf = (so.base.ope_1d(mult_sin_1d, 1, so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}

Val_domain Domain_polar_shell::div_sin_theta (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base = so.base ;	
	
	res.cf = (so.base.ope_1d(div_sin_1d, 1, so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}

Val_domain Domain_polar_shell::mult_r (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;

	res.cf = (so.base.ope_1d(mult_x_1d, 0, so.cf, res.base)*alpha + (so.cf)*beta) ;
	res.in_coef = true ;
	return res ;
}

Val_domain Domain_polar_shell::div_r (const Val_domain& so) const {
	Val_domain res (so / get_radius()) ;
	res.base = so.base ;
	return (res) ;
}

Val_domain Domain_polar_shell::laplacian (const Val_domain& so, int m) const {
  Val_domain derr (so.der_var(1)/alpha) ;
  Val_domain dert (so.der_var(2)) ;
  Val_domain res (derr.der_var(1)/alpha + div_r(2*derr + div_r(dert.der_var(2) + dert.mult_cos_theta().div_sin_theta()))) ;
  if (m!=0)
    res -= m * m * div_r(div_r(so.div_sin_theta().div_sin_theta())) ;
  return res ;
}

Val_domain Domain_polar_shell::laplacian2 (const Val_domain& so, int m) const {
  Val_domain derr (so.der_var(1)/alpha) ;
  Val_domain dert (so.der_var(2)) ;
  Val_domain res (derr.der_var(1)/alpha + div_r(derr + div_r(dert.der_var(2)))) ;
  if (m!=0)
    res -= m * m * div_r(div_r(so.div_sin_theta().div_sin_theta())) ;
  return res ;
}

Val_domain Domain_polar_shell::der_r (const Val_domain& so) const {
  return (so.der_var(1)/alpha) ;
}


Val_domain Domain_polar_shell::dt (const Val_domain& so) const {
  return (so.der_var(2)) ;
}


Val_domain Domain_polar_shell::div_xp1 (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;

	res.cf = (so.base.ope_1d(div_xp1_1d, 0, so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}

double Domain_polar_shell::integrale (const Val_domain& so) const {
  double res = 0 ;
  Val_domain integrant (mult_r(so)) ;
  integrant.get_coef() ;
  Array<double> cf (integrant.get_coef()) ;

  int baset = (*integrant.get_base().bases_1d[1]) (0) ;
  switch (baset) {
    case COS_ODD :
	break ;
    case SIN_EVEN :
	break ;
    case COS_EVEN : {
	// Only m=0 :
	double facttheta = M_PI ;
	int baser = (*integrant.get_base().bases_1d[0]) (0) ;
	switch (baser) {
	      case CHEB : {
		    for (int i=0 ; i<nbr_coefs(0) ; i+=2)
			res += facttheta * (1./double(i+1) - 1./double(i-1)) *cf(i,0) ;
		    break ;
		    }
	      case LEG : {
		   res += facttheta*2*cf(0,0) ;
		   break ;
		  }
	      default :
		  cerr << "Case not yet implemented in Domain_polar_shell::integrale" << endl ;
		  abort() ;
	}
	break ;
    }
    case SIN_ODD : {
	  for (int j=0 ; j<nbr_coefs(1) ; j++) {
	    double facttheta = 2./double(2*j+1) ;
	    int baser = (*integrant.get_base().bases_1d[0]) (0) ;
	  switch (baser) {
	      case CHEB : {
		    for (int i=0 ; i<nbr_coefs(0) ; i+=2)
			res += facttheta * (1./double(i+1) - 1./double(i-1)) *cf(i,j) ;
		    break ;
		    }
	      case LEG : {
		   res += facttheta*2*cf(0,j) ;
		   break ;
		  }
	      default :
		  cerr << "Case not yet implemented in Domain_polar_shell::integrale" << endl ;
		  abort() ;
	    }
      }
      break ;
    }
    default : 
      cerr << "Case not yet implemented in Domain_polar_shell::integrale" << endl ;
      abort() ;
  }
  // Phi contribution :
 // res *= 2*M_PI*alpha ;
  res *= alpha ;
  return res ;
}


double Domain_polar_shell::integ_volume (const Val_domain& so) const {
  return integrale(so) ;
}
}
