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
int div_sin_1d_pasprevu (Array<double>&) {
	cout << "div_sin_1d not implemented." << endl ;
	abort() ;
} 

int div_sin_1d_cos (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;

	int nr = tab.get_size(0) ;
	assert(nr%2==1) ;
	Array<double> res (nr) ;

	double somme = 0 ;
	res.set(nr-1) = 0. ;
	for (int k=nr-2; k>=0; k-=2) {
		somme -= 2*tab(k) ;
		res.set(k-1) = somme ;
		}
	somme = 0 ;
	for (int k=nr-3; k>=0; k-=2) {
		somme -= 2*tab(k+2) ;
		res.set(k+1) = somme ;
		}
	res.set(0) = 0 ;

	tab = res ;
	return SIN ;
}

int div_sin_1d_sin (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	assert (nr%2==1) ;
	Array<double> res (nr) ;

	double somme = 0 ;
	res.set(nr-1) = 0. ;
	for (int k=nr-2; k>=0; k-=2) {
		somme += 2*tab(k) ;
		res.set(k-1) = somme ;
		}
	somme = 0 ;
	for (int k=nr-3; k>=0; k-=2) {
		somme += 2*tab(k+2) ;
		res.set(k+1) = somme ;
		}
	res.set(0) *= 0.5 ;	

	tab = res ;
	return COS ;
}

int div_sin_1d_cos_even (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	Array<double> res (nr) ;

	double somme = 0 ;	
	res.set(nr-1) = 0. ;
	for (int k=nr-2; k>=0; k--) {
	    somme -= 2*tab(k+1) ;
	    res.set(k) = somme ;
	    }
	tab = res ;
	return SIN_ODD ;
}

int div_sin_1d_cos_odd (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	Array<double> res (nr) ;

	double somme = 0 ;	
	res.set(nr-1) = 0. ;
	for (int k=nr-2; k>=0; k--) {
	    somme -= 2*tab(k) ;
	    res.set(k) = somme ;
	}

	tab = res ;
	return SIN_EVEN ;
}


int div_sin_1d_sin_even (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	Array<double> res (nr) ;

	double somme = 0 ;	
	res.set(nr-1) = 0. ;
	for (int k=nr-2; k>=0; k--) {
	    somme += 2*tab(k+1) ;
	    res.set(k) = somme ;
	    }

	tab = res ;
	return COS_ODD ;
}

int div_sin_1d_sin_odd (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	Array<double> res (nr) ;

	double somme = 0 ;	
	res.set(nr-1) = 0. ;
	for (int k=nr-2; k>=0; k--) {
	    somme += 2*tab(k) ;
	    res.set(k) = somme ;
	    }
	res.set(0) *= 0.5 ;
	tab = res ;
	return COS_EVEN ;
}


int div_sin_1d (int base, Array<double>& tab) {
    static int (*div_sin_1d[NBR_MAX_BASE])(Array<double>&) ;
    static bool premier_appel = true ;
    
    // Premier appel
    if (premier_appel) {
	premier_appel = false ;

	for (int i=0; i<NBR_MAX_BASE; i++)
	    div_sin_1d[i] = div_sin_1d_pasprevu ;

	div_sin_1d[COS_EVEN] = div_sin_1d_cos_even ;
	div_sin_1d[COS_ODD] = div_sin_1d_cos_odd ;
	div_sin_1d[SIN_EVEN] = div_sin_1d_sin_even ;
	div_sin_1d[SIN_ODD] = div_sin_1d_sin_odd ;
	div_sin_1d[COS] = div_sin_1d_cos ;
	div_sin_1d[SIN] = div_sin_1d_sin ;
	}
	
        return div_sin_1d[base](tab) ;
}}
