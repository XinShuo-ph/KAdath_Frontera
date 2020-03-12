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

#include "array.hpp"
#include "adapted_polar.hpp"
#include "val_domain.hpp"
#include "array_math.hpp"
#include "scalar.hpp"
namespace Kadath {
int mult_cos_1d (int, Array<double>&) ;
int mult_sin_1d (int, Array<double>&) ;
int div_sin_1d (int, Array<double>&) ;
int div_cos_1d (int, Array<double>&) ;
int mult_x_1d (int, Array<double>&) ;

Val_domain Domain_polar_shell_outer_adapted::mult_cos_theta (const Val_domain& so) const {
	if (so.check_if_zero()) 
	  return so ;
	else {so.coef() ;
	Val_domain res(this) ;

	res.base = so.base ;	
	
	res.cf = new Array<double> (so.base.ope_1d(mult_cos_1d, 1, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;}
}

Val_domain Domain_polar_shell_outer_adapted::mult_sin_theta (const Val_domain& so) const {
	if (so.check_if_zero()) 
	  return so ;
	else {so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;
	
	res.cf = new Array<double> (so.base.ope_1d(mult_sin_1d, 1, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;}
}

Val_domain Domain_polar_shell_outer_adapted::div_sin_theta (const Val_domain& so) const {
	if (so.check_if_zero()) 
	  return so ;
	else {so.coef() ;
	Val_domain res(this) ;

	res.base = so.base ;	
	
	res.cf = new Array<double> (so.base.ope_1d(div_sin_1d, 1, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;}
}



Val_domain Domain_polar_shell_outer_adapted::der_r (const Val_domain& so) const {
  if (so.check_if_zero()) 
	  return so ;
	else {

	return (so.der_var(1) * 2. / (*outer_radius - inner_radius)) ;
	
	}
}

Val_domain Domain_polar_shell_outer_adapted::div_r (const Val_domain& so) const {
  if (so.check_if_zero()) 
	  return so ;
	else {

	return (so/ get_radius()) ;
	
	}
}

Val_domain Domain_polar_shell_outer_adapted::mult_r (const Val_domain& so) const {
  if (so.check_if_zero()) 
	  return so ;
	else {

	return (so* get_radius()) ;
	
	}
}

Val_domain Domain_polar_shell_outer_adapted::laplacian (const Val_domain& so, int m) const {
  Val_domain derr (so.der_r()) ;
  Val_domain dert (so.der_var(2)) ;
  Val_domain res (derr.der_r() + div_r(2*derr + div_r(dert.der_var(2) + dert.mult_cos_theta().div_sin_theta()))) ;
  if (m!=0)
    res -= m * m * div_r(div_r(so.div_sin_theta().div_sin_theta())) ;
  return res ;
}

Val_domain Domain_polar_shell_outer_adapted::laplacian2 (const Val_domain& so, int m) const {
  Val_domain derr (so.der_r()) ;
  Val_domain dert (so.der_var(2)) ;
  Val_domain res (derr.der_r() + div_r(derr + div_r(dert.der_var(2)))) ;
  if (m!=0)
    res -= m * m * div_r(div_r(so.div_sin_theta().div_sin_theta())) ;
  return res ;
}
}

