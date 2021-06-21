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
#include "val_domain.hpp"
#include "array.hpp"

namespace Kadath {
int mult_cos_1d (int, Array<double>&) ;
int mult_sin_1d (int, Array<double>&) ;
int div_sin_1d (int, Array<double>&) ;
int div_x_1d (int, Array<double>&) ;
int mult_x_1d (int, Array<double>&) ;
int div_xm1_1d (int, Array<double>&) ;

Val_domain Domain_polar_periodic_shell::mult_cos_theta (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;
	res.base = so.base ;
	res.cf = new Array<double> (so.base.ope_1d(mult_cos_1d, 1, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}

Val_domain Domain_polar_periodic_shell::mult_sin_theta (const Val_domain& so) const {

	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;
	res.cf = new Array<double> (so.base.ope_1d(mult_sin_1d, 1, *so.cf, res.base)) ;
	res.in_coef = true ;	
	return res ;
}

Val_domain Domain_polar_periodic_shell::div_sin_theta (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base = so.base ;	
	
	res.cf = new Array<double> (so.base.ope_1d(div_sin_1d, 1, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}


Val_domain Domain_polar_periodic_shell::mult_r (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;
	
	res.cf = new Array<double> (so.base.ope_1d(mult_x_1d, 0, *so.cf, res.base)*alpha + (*so.cf)*beta) ;
	res.in_coef = true ;
	return res ;
}


Val_domain Domain_polar_periodic_shell::div_r (const Val_domain& so) const {
	Val_domain res (so / get_radius()) ;
	res.base = so.base ;
	return (res) ;
}


Val_domain Domain_polar_periodic_shell::div_1mrsL (const Val_domain& so) const {
   if (so.check_if_zero())
      return so;
	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;

	res.cf = new Array<double> (-(alpha+beta)/alpha*so.base.ope_1d(div_xm1_1d, 0, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}



Val_domain Domain_polar_periodic_shell::der_r (const Val_domain& so) const {
  return (so.der_var(1)/alpha) ;
}

Val_domain Domain_polar_periodic_shell::dt (const Val_domain& so) const {
  return (so.der_var(2)) ;
}

Val_domain Domain_polar_periodic_shell::dtime (const Val_domain& so) const {
  return (so.der_var(3)*ome) ;
}


Val_domain Domain_polar_periodic_shell::mult_cos_time (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;
	res.base = so.base ;
	res.cf = new Array<double> (so.base.ope_1d(mult_cos_1d, 2, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}

Val_domain Domain_polar_periodic_shell::mult_sin_time (const Val_domain& so) const {

	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;
	res.cf = new Array<double> (so.base.ope_1d(mult_sin_1d, 2, *so.cf, res.base)) ;
	res.in_coef = true ;	
	return res ;
}
}

