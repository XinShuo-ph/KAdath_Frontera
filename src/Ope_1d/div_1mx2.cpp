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

#include "base_spectral.hpp"
#include "headcpp.hpp"
#include "matrice.hpp"
#include "array.hpp"
namespace Kadath {
int div_1mx2_1d_pasprevu (Array<double>&) {
	cout << "div_1mx2_1d not implemented." << endl ;
	abort() ;
} 


int div_1mx2_1d_cheb (Array<double>& so) {
	assert (so.get_ndim()==1) ;
	int nr = so.get_size(0) ;
	Array<double> res (nr) ;
	res = 0 ;
	for (int n=0 ; n<nr ; n++)
		for (int p=n+2 ; p<nr ; p+=2)
			res.set(n) += -2*(p-n)*so(p) ;
	res.set(0) *= 0.5 ;
	
	so = res ;
	return CHEB ;
}

int div_1mx2_1d_cheb_even (Array<double>& so) {
	assert (so.get_ndim()==1) ;
	int nr = so.get_size(0) ;
	Array<double> res (nr) ;
	res = 0 ;
	for (int n=0 ; n<nr ; n++)
		for (int p=n+1 ; p<nr ; p++)
			res.set(n) += -4*(p-n)*so(p) ;
	res.set(0) *= 0.5 ;
	
	so = res ;
	return CHEB_EVEN ;
}

int div_1mx2_1d_cheb_odd (Array<double>& so) {
	assert (so.get_ndim()==1) ;
	int nr = so.get_size(0) ;
	Array<double> res (nr) ;
	res = 0 ;
	for (int n=0 ; n<nr ; n++)
		for (int p=n+1 ; p<nr ; p++)
			res.set(n) += -4*(p-n)*so(p) ;
	
	so = res ;
	return CHEB_ODD ;
}

int div_1mx2_1d (int base, Array<double>& so) {
    static int (*div_1mx2_1d[NBR_MAX_BASE])(Array<double>&) ;
    static bool premier_appel = true ;
    
    // Premier appel
    if (premier_appel) {
	premier_appel = false ;

	for (int i=0; i<NBR_MAX_BASE; i++)
	    div_1mx2_1d[i] = div_1mx2_1d_pasprevu ;

	div_1mx2_1d[CHEB] = div_1mx2_1d_cheb ;
	div_1mx2_1d[CHEB_EVEN] = div_1mx2_1d_cheb_even ;
	div_1mx2_1d[CHEB_ODD] = div_1mx2_1d_cheb_odd ;
	}
	
        return div_1mx2_1d[base](so) ;
}}
