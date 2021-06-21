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
#include "array.hpp"
namespace Kadath {
int mult_sin_1d_pasprevu (Array<double>&) {
	cout << "mult_sin_1d not implemented." << endl ;
	abort() ;
} 

int mult_sin_1d_cossin (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	assert (nr%2==0) ;
	Array<double> res (nr) ;

	res.set(0) = 0.5*tab(3) ;
	res.set(1) = 0 ;
	res.set(2) = 0.5*tab(5) ;
	res.set(3) = tab(0) - 0.5*tab(4) ;
	for (int k=4; k<nr-2; k+=2) {
            // Cosines :
	    res.set(k) = 0.5*(tab(3+k) - tab(k-1)) ;
	    // Sines
	    res.set(k+1) = 0.5*(tab(k-2) - tab(k+2)) ;
	    }
	
	res.set(nr-2) = -0.5*tab(nr-3) ;
	res.set(nr-1) = 0. ;

	tab = res ;
	return COSSIN ;
}

int mult_sin_1d_cos (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	Array<double> res (nr) ;

	res.set(0) = 0. ;
	res.set(1) = tab(0) - 0.5*tab(2) ;
	for (int k=2; k<nr-1; k++)
	    // Sines
	    res.set(k) = 0.5*(tab(k-1) - tab(k+1)) ;
	//res.set(nr-1) = 0.5*tab(nr-2) ;
	res.set(nr-1) = 0 ;

	tab = res ;
	return SIN ;
}

int mult_sin_1d_sin (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	Array<double> res (nr) ;

	res.set(0) = 0.5*tab(1) ;
	res.set(1) = 0.5*tab(2) ;
	for (int k=2; k<nr-1; k++)
	    // Sines
	    res.set(k) = 0.5*(-tab(k-1) + tab(k+1)) ;
	//res.set(nr-1) = -0.5*tab(nr-2) ;
	res.set(nr-1) = 0 ;

	tab = res ;
	return COS ;
}

int mult_sin_1d_cos_even (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	Array<double> res (nr) ;

	res.set(0) = tab(0) - 0.5*tab(1) ;
	for (int k=1; k<nr-1; k++)
	    // Sines
	    res.set(k) = 0.5*(tab(k) - tab(k+1)) ;
	//res.set(nr-1) = 0.5*tab(nr-1) ;
	res.set(nr-1) = 0 ;

	tab = res ;
	return SIN_ODD ;
}

int mult_sin_1d_cos_odd (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	Array<double> res (nr) ;

	res.set(0) = 0. ;
	for (int k=1; k<nr; k++)
	    // Sines
	    res.set(k) = 0.5*(tab(k-1) - tab(k)) ;

	tab = res ;
	return SIN_EVEN ;
}

int mult_sin_1d_sin_even (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	Array<double> res (nr) ;

	res.set(0) = 0.5*tab(1) ;
	for (int k=1; k<nr-1; k++)
	    // Sines
	    res.set(k) = 0.5*(-tab(k) + tab(k+1)) ;
	//res.set(nr-1) = -0.5*tab(nr-1) ;
	res.set(nr-1) = 0 ;

	tab = res ;
	return COS_ODD ;
}

int mult_sin_1d_sin_odd (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	Array<double> res (nr) ;

	res.set(0) = 0.5*tab(0) ;
	for (int k=1; k<nr ; k++)
	    // Sines
	    res.set(k) = 0.5*(-tab(k-1) + tab(k)) ;

	tab = res ;
	return COS_EVEN ;
}

int mult_sin_1d (int base, Array<double>& tab) {
    static int (*mult_sin_1d[NBR_MAX_BASE])(Array<double>&) ;
    static bool premier_appel = true ;
    
    // Premier appel
    if (premier_appel) {
	premier_appel = false ;

	for (int i=0; i<NBR_MAX_BASE; i++)
	    mult_sin_1d[i] = mult_sin_1d_pasprevu ;

	mult_sin_1d[COSSIN] = mult_sin_1d_cossin ;
	mult_sin_1d[COS_EVEN] = mult_sin_1d_cos_even ;
	mult_sin_1d[COS_ODD] = mult_sin_1d_cos_odd ;
	mult_sin_1d[SIN_EVEN] = mult_sin_1d_sin_even ;
	mult_sin_1d[SIN_ODD] = mult_sin_1d_sin_odd ;
	mult_sin_1d[COS] = mult_sin_1d_cos ;
	mult_sin_1d[SIN] = mult_sin_1d_sin ;
	}
	
        return mult_sin_1d[base](tab) ;
}}
