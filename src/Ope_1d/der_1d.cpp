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
#include "array_math.hpp"
namespace Kadath {
int der_1d_pasprevu (Array<double>&) {
	cout << "Der_1d not implemented." << endl ;
	abort() ;
} 

	// CHEB :
int der_1d_cheb (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	
	Array<double> res (nr) ;
	res = 0 ;
	for (int i=0 ;i<nr ; i++)
	    for (int p=i+1 ; p<nr ; p+=2)
	    	res.set(i) += p*tab(p) ;
	res *= 2 ;
	res.set(0) /= 2 ;
	tab = res ;
	
	
	// Output base :
	return CHEB ;
}
	// CHEB EVEN:
int der_1d_cheb_even (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;
	res = 0 ;
	for (int i=0 ;i<nr ; i++)
	    for (int p=i+1 ; p<nr ; p++)
	    	res.set(i) += 2*p*tab(p) ;
	res.set(nr-1) = 0 ;
	res *= 2 ;
	tab = res ;
	
	// Output base :
	return CHEB_ODD ;
}

	// CHEB ODD:
int der_1d_cheb_odd (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;
	res = 0 ;
	for (int i=0 ;i<nr ; i++)
	    for (int p=i ; p<nr ; p++)
	    	res.set(i) += (2*p+1)*tab(p) ;
	res *= 2 ;
	res.set(0) /= 2 ;
	tab = res ;
	
	// Output base :
	return CHEB_EVEN ;
}

	// LEG :
int der_1d_leg (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;
	res = 0 ;
	for (int i=0 ;i<nr ; i++)
	    for (int p=i+1 ; p<nr ; p+=2)
	    	res.set(i) += (2*i+1)*tab(p) ;
	tab = res ;
	
	// Output base :
	return LEG ;
}
	// LEG EVEN:
int der_1d_leg_even (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;
	res = 0 ;
	for (int i=0 ;i<nr ; i++)
	    for (int p=i+1 ; p<nr ; p++)
	    	res.set(i) += (4*i+3)*tab(p) ;
	res.set(nr-1) = 0 ;
	tab = res ;	

	// Output base :
	return LEG_ODD ;
}

	// LEG ODD:
int der_1d_leg_odd (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;
	res = 0 ;
	for (int i=0 ;i<nr ; i++)
	    for (int p=i ; p<nr ; p++)
	    	res.set(i) += (4*i+1)*tab(p) ;
	tab = res ;
	
	// Output base :
	return LEG_EVEN ;
}

	// COSSIN:
int der_1d_cossin (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;
	// cosines
	for (int i=0 ;i<nr-1 ; i+=2)
	    res.set(i+1) = - (i/2)*tab(i) ;
	// Sines
	for (int i=1 ;i<nr ; i+=2)
	    res.set(i-1) = ((i-1)/2)*tab(i) ;
	res.set(nr-1) = 0 ;
	tab = res ;
	
	// Output base :
	return COSSIN ;
}

	// COSSIN EVEN :
int der_1d_cossin_even (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;
	// cosines
	for (int i=0 ;i<nr-1 ; i+=2)
	    res.set(i+1) = - i*tab(i) ;
	// Sines
	for (int i=1 ;i<nr ; i+=2)
	    res.set(i-1) = (i-1)*tab(i) ;
	tab = res ;
	
	// Output base :
	return COSSIN_EVEN ;
}

	// COSSIN ODD :
int der_1d_cossin_odd (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;
	// cosines
	for (int i=0 ;i<nr-1 ; i+=2)
	    res.set(i+1) = - (i+1)*tab(i) ;
	// Sines
	for (int i=1 ;i<nr ; i+=2)
	    res.set(i-1) = i*tab(i) ;
	tab = res ;
	
	// Output base :
	return COSSIN_ODD ;
}

	// COS:
int der_1d_cos (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;

	for (int i=0 ;i<nr ; i++)
	    res.set(i) = - i*tab(i) ;
	tab = res ;
	
	// Output base :
	return SIN ;
}


	// SIN:
int der_1d_sin(Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;

	for (int i=0 ;i<nr ; i++)
	    res.set(i) = i*tab(i) ;
	tab = res ;
	
	// Output base :
	return COS ;
}



	// COS EVEN:
int der_1d_cos_even (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;

	for (int i=0 ;i<nr ; i++)
	    res.set(i) = -2*i*tab(i) ;
	res.set(nr-1) = 0 ;
	tab = res ;
	
	// Output base :
	return SIN_EVEN ;
}

	
	// COS ODD:
int der_1d_cos_odd (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;

	for (int i=0 ;i<nr ; i++)
	    res.set(i) = -(2*i+1)*tab(i) ;
	tab = res ;
	
	// Output base :
	return SIN_ODD ;
}

	


	// SIN EVEN:
int der_1d_sin_even (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;

	for (int i=0 ;i<nr ; i++)
	    res.set(i) = 2*i*tab(i) ;
	tab = res ;
	
	// Output base :
	return COS_EVEN ;
}

	
	// SIN ODD:
int der_1d_sin_odd (Array<double>& tab) {
	
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	Array<double> res (nr) ;

	for (int i=0 ;i<nr ; i++)
	    res.set(i) = (2*i+1)*tab(i) ;
	res.set(nr-1) = 0 ;
	tab = res ;
	
	// Output base :
	return COS_ODD ;
}

	

int der_1d (int base, Array<double>& tab) {
    static int (*der_1d[NBR_MAX_BASE])(Array<double>&) ;
    static bool premier_appel = true ;
    
    // Premier appel
    if (premier_appel) {
	premier_appel = false ;

	for (int i=0; i<NBR_MAX_BASE; i++)
	    der_1d[i] = der_1d_pasprevu ;

	der_1d[CHEB] = der_1d_cheb ;
	der_1d[CHEB_EVEN] = der_1d_cheb_even ;
	der_1d[CHEB_ODD] = der_1d_cheb_odd ;
	der_1d[COSSIN] = der_1d_cossin ;
	der_1d[COS_EVEN] = der_1d_cos_even ;
	der_1d[COS_ODD] = der_1d_cos_odd ;
	der_1d[SIN_EVEN] = der_1d_sin_even ;
	der_1d[SIN_ODD] = der_1d_sin_odd ;
	der_1d[COS] = der_1d_cos ;
	der_1d[SIN] = der_1d_sin ;
	der_1d[LEG] = der_1d_leg ;
	der_1d[LEG_EVEN] = der_1d_leg_even ;
	der_1d[LEG_ODD] = der_1d_leg_odd;
	der_1d[COSSIN_EVEN] = der_1d_cossin_even ;
	der_1d[COSSIN_ODD] = der_1d_cossin_odd ;
	}
	
        return der_1d[base](tab) ;
}}
