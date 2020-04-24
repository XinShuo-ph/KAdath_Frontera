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
#include "val_domain.hpp"
#include "array_math.hpp"
namespace Kadath {
int mult_x_1d (int, Array<double>&) ;

Val_domain Domain_spheric_time_shell::mult_r (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;

	res.cf = new Array<double> (so.base.ope_1d(mult_x_1d, 0, *so.cf, res.base)*alpha + (*so.cf)*beta) ;
	res.in_coef = true ;
	return res ;
}

Val_domain Domain_spheric_time_shell::div_r (const Val_domain& so) const {
	Val_domain res (so / get_radius()) ;
	res.base = so.base ;
	return (res) ;
}

Val_domain Domain_spheric_time_shell::der_r (const Val_domain& so) const {
  return (so.der_var(1)/alpha) ;
}

Val_domain Domain_spheric_time_shell::ddt (const Val_domain& so) const {
  return (so.der_var(2).der_var(2)*4./(tmax-tmin)/(tmax-tmin)) ;
}

Val_domain Domain_spheric_time_shell::dt (const Val_domain& so) const {
  return (so.der_var(2)*2./(tmax-tmin)) ;
}

Val_domain Domain_spheric_time_shell::laplacian (const Val_domain& so, int m) const {

#ifndef REMOVE_ALL_CHECKS
  if (m!=0) {
      cerr << "Laplacian only defiend for m=0 for Domain_spheric_shell" << endl ;
      abort() ;
  }
#endif
  Val_domain dr (so.der_var(1)/alpha) ;

  Val_domain res (dr.der_var(1)/alpha + div_r(2*dr)) ;
  return res ;
}}
