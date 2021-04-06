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
#include "val_domain.hpp"
#include "array_math.cpp"
#include <gsl/gsl_sf_bessel.h>

namespace Kadath {

void Val_domain::operator+= (const Val_domain& so) {
	*this = *this + so ;
}

void Val_domain::operator-= (const Val_domain& so) {
	*this = *this - so ;
}

void Val_domain::operator*= (const Val_domain& so) {
	*this = *this * so ;
}

void Val_domain::operator/= (const Val_domain& so) {
	*this = *this / so ;
}

void Val_domain::operator+= (double xx) {
	*this = *this + xx ;
}

void Val_domain::operator-= (double xx) {
	*this = *this - xx ;
}
	

void Val_domain::operator*= (double xx) {
	*this = *this * xx ;
}
	

void Val_domain::operator/= (double xx) {
	*this = *this / xx ;
}
		
	
Val_domain sin(const Val_domain& so) {
	if (so.is_zero)
		return so ;
	else {
		so.coef_i() ;
		Val_domain res (so.zone) ;
		res.set_in_conf() ;
		res.c = new Array<double> (sin(*so.c)) ;
		//res.std_base();
		return res ;
	}
}


Val_domain cos(const Val_domain& so ) {

	if (so.is_zero) {
		Val_domain res(so) ;
		res.is_zero = false ;
		res.allocate_conf() ;
		*res.c = 1. ;
		//res.std_base();
		return res ;
	}
	else {
		so.coef_i() ;	
		Val_domain res (so.zone) ;
		res.set_in_conf() ;
		res.c = new Array<double> (cos(*so.c)) ;
		//res.std_base();
		return res ;
	}
}

Val_domain cosh(const Val_domain& so ) {

	if (so.is_zero) {
		Val_domain res(so) ;
		res.is_zero = false ;
		res.allocate_conf() ;
		*res.c = 1. ;
		//res.std_base();
		return res ;
	}
	else {
		so.coef_i() ;	
		Val_domain res (so.zone) ;
		res.set_in_conf() ;
		res.c = new Array<double> (cosh(*so.c)) ;
		//res.std_base();
		return res ;
	}
}

Val_domain sinh(const Val_domain& so ) {

	if (so.is_zero) {
		Val_domain res(so) ;
		res.is_zero = false ;
		res.allocate_conf() ;
		*res.c = 1. ;
		//res.std_base();
		return res ;
	}
	else {
		so.coef_i() ;	
		Val_domain res (so.zone) ;
		res.set_in_conf() ;
		res.c = new Array<double> (sinh(*so.c)) ;
		//res.std_base();
		return res ;
	}
}

Val_domain operator+ (const Val_domain& so) {
	Val_domain res (so) ;
	return res ;
}

Val_domain operator- (const Val_domain& so) {
	Val_domain res (so) ;
	if (so.in_conf)
		*res.c *= -1 ;
	if (so.in_coef) 
		*res.cf *= -1 ;
	return res ;
}


Val_domain operator+ (const Val_domain& a, const Val_domain& b) {
	assert (a.zone == b.zone) ;

	if (a.is_zero)
		return b ;
	else if (b.is_zero)
			return a ;
	else {
	Val_domain res (a.zone) ;

	if (a.in_conf) {
		if (!b.in_conf)
			b.coef_i() ;
		res.c = new Array<double> (*a.c + *b.c) ;
		res.in_conf=true ;
		}
	else if (b.in_conf) {
		if (!a.in_conf)
			a.coef_i() ;
		res.c = new Array<double> (*a.c + *b.c) ;
		res.in_conf=true ;
		}
	
	if (!res.in_conf) {
	if (a.in_coef) {
		if (!b.in_coef)
			b.coef() ;
		res.cf = new Array<double> (*a.cf + *b.cf) ;
		res.in_coef=true ;
		}
	else if (b.in_coef) {
		if (!a.in_coef)
			a.coef() ;
		res.cf = new Array<double> (*a.cf + *b.cf) ;
		res.in_coef=true ;
		}	
	}
	if ((a.base.is_def()) && (b.base.is_def())) {
		assert (a.base==b.base) ;
		res.base = a.base ;
	}

	return res ;
	}
}


Val_domain operator+ (const Val_domain& so, double x) {
	if (so.is_zero) {
		Val_domain res (so.zone) ;
		res.allocate_conf() ;
		*res.c = x ;
		res.std_base() ;
		return res ;
	}
	else {
		so.coef_i() ;
		Val_domain res (so.zone) ;
		res.set_in_conf() ;
		res.c = new Array<double> (*so.c + x) ;
		res.base = so.base ;
		return res ;
	}
}


Val_domain operator+ (double x, const Val_domain& so) {
	return so+x ;
}

Val_domain operator- (const Val_domain& a, const Val_domain& b) {
	assert (a.zone == b.zone) ;

	if (a.is_zero)
		return -b ;
	else if (b.is_zero)
		return a ;
	else {
	Val_domain res (a.zone) ;
	if (a.in_conf) {
		if (!b.in_conf)
			b.coef_i() ;
		res.c = new Array<double> (*a.c - *b.c) ;
		res.in_conf=true ;
		}
	else if (b.in_conf) {
		if (!a.in_conf)
			a.coef_i() ;
		res.c = new Array<double> (*a.c - *b.c) ;
		res.in_conf=true ;
		}
	if (!res.in_conf) {
	if (a.in_coef) {
		if (!b.in_coef)
			b.coef() ;
		assert (a.base == b.base) ;
		res.cf = new Array<double> (*a.cf - *b.cf) ;
		res.in_coef=true ;
		}
	else if (b.in_coef) {
		if (!a.in_coef)
			a.coef() ;
		res.cf = new Array<double> (*a.cf - *b.cf) ;
		res.in_coef=true ;
		}
	}
	if ((a.base.is_def()) && (b.base.is_def())) {
		assert (a.base==b.base) ;
		res.base = a.base ;
	}

	return res ;
	}
}

Val_domain operator- (const Val_domain& so, double x) {	

	if (so.is_zero) {
		Val_domain res (so.zone) ;
		res.allocate_conf() ;
		*res.c = -x ;
		res.std_base() ;
		return res ;
	}
	else {
		so.coef_i() ;
		Val_domain res (so.zone) ;
		res.set_in_conf() ;
		res.c = new Array<double> (*so.c - x) ;
		res.base = so.base ;
		return res ;
	}
}

Val_domain operator- (double x, const Val_domain& so) {
	if (so.is_zero) {
		Val_domain res (so.zone) ;
		res.allocate_conf() ;
		*res.c = x ;
		res.std_base() ;
		return res ;
	}
	else {
		so.coef_i() ;
		Val_domain res (so.zone) ;
		res.set_in_conf() ;
		res.c = new Array<double> (x-*so.c) ;
		res.base = so.base ;
		return res ;
	}
}

Val_domain operator* (const Val_domain& a, const Val_domain& b) {
	if (a.is_zero)
		return a ;
	else if (b.is_zero)
		return b ;
		else {
	a.coef_i() ;
	b.coef_i() ;
        assert (a.zone==b.zone) ;
	Val_domain res (a.zone) ;
	res.allocate_conf() ;
	*res.c = (*a.c)* (*b.c) ;
	res.base = a.zone->mult(a.base, b.base) ;
	return res ;
	}
}
	

Val_domain operator* (const Val_domain& so, double x) {
	if (so.is_zero)
		return so ;
	else {
	Val_domain res (so) ;
	if (res.in_conf)
		*res.c *= x ;
	if (res.in_coef)
		*res.cf *= x ;
	return res ;
	}
}

Val_domain operator* (double x, const Val_domain& so) {
	return (so*x);
}

Val_domain operator* (const Val_domain& so, int nn) {
	if (so.is_zero)
		return so ;
	else 
		if (nn==0) {
			Val_domain res(so, false) ;
			res.set_zero() ;
			return res ;
		} 
		else {
		Val_domain res (so) ;
		if (res.in_conf)
			*res.c *= nn ;
		if (res.in_coef)
			*res.cf *= nn ;
		return res ;
		}
}

Val_domain operator* (int nn, const Val_domain& so) {
	return (so*nn);
}

Val_domain operator* (const Val_domain& so, long int nn) {
	if (so.is_zero)
		return so ;
	else 
		if (nn==0) {
			Val_domain res(so, false) ;
			res.set_zero() ;
			return res ;
		} 
		else {
		Val_domain res (so) ;
		if (res.in_conf)
			*res.c *= static_cast<double>(nn) ;
		if (res.in_coef)
			*res.cf *= static_cast<double>(nn) ;
		return res ;
		}
}

Val_domain operator* (long int nn, const Val_domain& so) {
	return (so*nn);
}
Val_domain operator/ (const Val_domain& a, const Val_domain& b) {
	if (a.is_zero)
		return a ;
	else if (b.is_zero) {
		cerr << "Division by zero" << endl ;
		abort() ;
		}
	else {
	a.coef_i() ;
	b.coef_i() ;
        assert (a.zone==b.zone) ;
	Val_domain res (a.zone) ;
	res.allocate_conf() ;
	*res.c = (*a.c)/ (*b.c) ;
	res.base = a.zone->mult(a.base, b.base) ;
	return res ;
	}
}

Val_domain operator/ (const Val_domain& so, double x) {
	if (so.is_zero)
		return so ;
	else {
	Val_domain res (so) ;
	if (res.in_conf)
		*res.c /= x ;
	if (res.in_coef)
		*res.cf /= x ;
	return res ;
	}
}

Val_domain operator/ (double x, const Val_domain& so) {
	if (so.is_zero) {
		cerr << "Division by zero" << endl ;
		abort() ;
	}
	else {
	so.coef_i() ;
	Val_domain res (so.zone) ;
	res.set_in_conf() ;
	res.c = new Array<double> (x/ *so.c) ;
	// No sense if the base is not the standard one !
	res.std_base() ;
	return res ;
	}
}

Val_domain pow (const Val_domain& so, int n) {
	if (so.is_zero)
		return so ;
	else {
	if (n<=0) {
		so.coef_i() ;
		Val_domain res (so.zone) ;
		res.set_in_conf() ;
		res.c = new Array<double> (pow(*so.c,n)) ;
		return res ;
	}
	else {
		Val_domain res (so) ;
		for (int i=1 ; i<n ; i++)
			res *= so ;
		return res ;
		}
	}
}

Val_domain pow (const Val_domain& so, double nn) {
	if (so.is_zero) 
		return so ;
	else {
	so.coef_i() ;
	Val_domain res (so.zone) ;
	res.set_in_conf() ;
	res.c = new Array<double> (pow(*so.c,nn)) ;
	return res ;
	}
}

Val_domain sqrt (const Val_domain& so) {
	if (so.is_zero) 
		return so ;
	else {
	so.coef_i() ;
	Val_domain res (so.zone) ;
	res.set_in_conf() ;
	res.c = new Array<double> (sqrt(*so.c)) ;
	// Provisory ? same base as so...
	res.base = so.base ;
	return res ;
	}
}

Val_domain exp (const Val_domain& so) {
	if (so.is_zero) {
		Val_domain res(so) ;
		res.is_zero = false ;
		res.allocate_conf() ;
		*res.c = 1. ;
		return res ;
	}
	else {
	so.coef_i() ;
	Val_domain res (so.zone) ;
	res.set_in_conf() ;
	res.c = new Array<double> (exp(*so.c)) ;
	// Provisory ? same base as so...
	res.base = so.base ;
	return res ;
	}
}

Val_domain log (const Val_domain& so) {
	
	so.coef_i() ;
	Val_domain res (so.zone) ;
	res.set_in_conf() ;
	res.c = new Array<double> (log(*so.c)) ;
	// Provisory ? same base as so...
	res.base = so.base ;
	return res ;
}

Val_domain atanh (const Val_domain& so) {
	
	so.coef_i() ;
	Val_domain res (so.zone) ;
	res.set_in_conf() ;
	res.c = new Array<double> (atanh(*so.c)) ;
	// Provisory ? same base as so...
	res.base = so.base ;
	return res ;
}

double diffmax (const Val_domain& aa, const Val_domain& bb) {
	aa.coef_i() ;
	bb.coef_i() ;
	return diffmax (*aa.c, *bb.c) ;
}

Val_domain bessel_jl (const Val_domain& so, int l) {
  
    Val_domain res (so.get_domain()) ;
    res.allocate_conf() ;
    Index pos (res.get_domain()->get_nbr_points()) ;
    do  {
      res.set(pos) = gsl_sf_bessel_jl (l, so(pos)) ;
    }
    while (pos.inc()) ;
    res.std_base() ;
    return res ;
}


Val_domain bessel_yl (const Val_domain& so, int l) {
  
    Val_domain res (so.get_domain()) ;
    res.allocate_conf() ;
    Index pos (res.get_domain()->get_nbr_points()) ;
    do  {
      res.set(pos) = gsl_sf_bessel_yl (l, so(pos)) ;
    }
    while (pos.inc()) ;
    res.std_base() ;
    return res ;
}

Val_domain bessel_djl (const Val_domain& so, int l) {
  
    Val_domain res (so.get_domain()) ;
    res.allocate_conf() ;
    Index pos (res.get_domain()->get_nbr_points()) ;
    do  {  
      if (l==0)
	res.set(pos) = cos(so(pos))/so(pos) - sin(so(pos))/so(pos)/so(pos) ;
      else
      res.set(pos) = gsl_sf_bessel_jl (l-1, so(pos)) - (l+1) * gsl_sf_bessel_jl (l, so(pos))/so(pos) ;
    }
    while (pos.inc()) ;
    res.std_base() ;
    return res ;
}


Val_domain bessel_dyl (const Val_domain& so, int l) {
  
    Val_domain res (so.get_domain()) ;
    res.allocate_conf() ;
    Index pos (res.get_domain()->get_nbr_points()) ;
    do  {
      if (l==0)
	res.set(pos) = sin(so(pos))/so(pos) + cos(so(pos))/so(pos)/so(pos) ;
      else
      res.set(pos) = gsl_sf_bessel_yl (l-1, so(pos)) - (l+1) * gsl_sf_bessel_yl (l, so(pos))/so(pos) ;
    }
    while (pos.inc()) ;
    res.std_base() ;
    return res ;
}

Val_domain atan (const Val_domain& so) {
	
	so.coef_i() ;
	Val_domain res (so.zone) ;
	res.set_in_conf() ;
	res.c = new Array<double> (atan(*so.c)) ;
	// Provisory ? same base as so...
	res.base = so.base ;
	return res ;
}

double maxval (const Val_domain& target)
{
   Val_domain so(target);
   double res(0);
   double value(0.0);
   if (!so.check_if_zero()) 
   {
      so.coef_i() ;
      Index pos(so.get_domain()->get_nbr_points()) ;
      do 
      {
         value = so(pos);
         if (fabs(value) > res) res = fabs(value);
      }while (pos.inc());
   }
   return res ;
}
}

