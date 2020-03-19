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

#include "point.hpp"
#include "index.hpp"
#include "base_spectral.hpp"
#include "array.hpp"
namespace Kadath {
double summation_1d_pasprevu (double, const Array<double>&) {
	cout << "Summation_1d not implemented." << endl ;
	abort() ;
} 

double summation_1d_cheb (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	double tm2 = 1. ;
	double tm1 = xx ;
	double tm ;
	
	double res = tab(0)*tm2 ;
	if (nr>1)
		res += tab(1)*tm1 ;
	for (int i=2 ; i<nr ; i++) {
		tm = 2*xx*tm1 - tm2 ;
		res += tab(i)*tm ;
		tm2 = tm1 ;
		tm1 = tm ;
	}
	return res ;
}

double summation_1d_cheb_even (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	double tm2 = 1. ;
	double tm1 = xx ;
	double tm ;
	
	double res = tab(0)*tm2 ;
	for (int i=1 ; i<nr ; i++) {
		tm = 2*xx*tm1 - tm2 ;
		res += tab(i)*tm ;
		tm2 = tm1 ;
		tm1 = tm ;
		tm = 2*xx*tm1 - tm2 ;
		tm2 = tm1 ;
		tm1 = tm ;
	}
	return res ;
}

double summation_1d_cheb_odd (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	double tm2 = 1. ;
	double tm1 = xx ;
	double tm ;
	
	
	double res = tab(0)*tm1 ;
	for (int i=1 ; i<nr ; i++) {
		tm = 2*xx*tm1 - tm2 ;
		tm2 = tm1 ;
		tm1 = tm ;
		tm = 2*xx*tm1 - tm2 ;
		res += tab(i)*tm ;
		tm2 = tm1 ;
		tm1 = tm ;
	}
	return res ;
}

double summation_1d_leg(double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	double plm2 = 1. ;
	double plm1 = xx ;
	double pl ;
	
	double res = tab(0)*plm2 ;
	if (nr>1)
		res += tab(1)*plm1 ;
	for (int l=2 ; l<nr ; l++) {
		pl = ((2*l-1)*xx*plm1 - (l-1)*plm2)/l ;
		res += tab(l)*pl ;
		plm2 = plm1 ;
		plm1 = pl ;
	}
	return res ;
}

double summation_1d_leg_even (double xx, const Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	double plm2 = 1. ;
	double plm1 = xx ;
	double pl ;
	
	double res = tab(0)*plm2 ;
	int courant ;
	for (int l=1 ; l<nr ; l++) {
		courant = 2*l ;
		pl = ((2*courant-1)*xx*plm1 - (courant-1)*plm2)/courant ;
		res += tab(l)*pl ;
		plm2 = plm1 ;
		plm1 = pl ;
		courant ++ ;
		pl = ((2*courant-1)*xx*plm1 - (courant-1)*plm2)/courant ;
		plm2 = plm1 ;
		plm1 = pl ;
	}

	return res ;
}


double summation_1d_leg_odd (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;
	double plm2 = 1. ;
	double plm1 = xx ;
	double pl ;
	
	double res = tab(0)*plm1 ;
	int courant ;
	for (int l=1 ; l<nr ; l++) {
		courant = 2*l ;
		pl = ((2*courant-1)*xx*plm1 - (courant-1)*plm2)/courant ;
		plm2 = plm1 ;
		plm1 = pl ;
		courant ++ ;
		pl = ((2*courant-1)*xx*plm1 - (courant-1)*plm2)/courant ;
		res += tab(l)*pl ;
		plm2 = plm1 ;
		plm1 = pl ;
	}
	return res ;
}

double summation_1d_cossin (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	double res = 0 ;
	for (int i=0 ; i<nbr ; i++)
		res += (i%2==0) ? tab(i)*cos((i/2)*xx) : tab(i)*sin((i-1)/2*xx) ;
	return res ;
}


double summation_1d_cossin_even (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	double res = 0 ;
	for (int i=0 ; i<nbr ; i++) 
		res += (i%2==0) ? tab(i)*cos(i*xx) : tab(i)*sin((i-1)*xx);
	return res ;
}

double summation_1d_cossin_odd (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	double res = 0 ;
	for (int i=0 ; i<nbr ; i++) 
		res += (i%2==0) ? tab(i)*cos((i+1)*xx) : tab(i)*sin(i*xx);
	return res ;
}

double summation_1d_cos (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	double res = 0 ;
	for (int i=0 ; i<nbr ; i++) 
		res += tab(i)*cos(i*xx) ;
	return res ;
}

double summation_1d_sin (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	double res = 0 ;
	for (int i=1 ; i<nbr ; i++) 
		res += tab(i)*sin(i*xx) ;
	return res ;
}

double summation_1d_cos_even (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	double res = 0 ;
	for (int i=0 ; i<nbr ; i++) 
		res += tab(i)*cos(2*i*xx) ;
	return res ;
}

double summation_1d_cos_odd (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	double res = 0 ;
	for (int i=0 ; i<nbr ; i++) 
		res += tab(i)*cos((2*i+1)*xx) ;
	return res ;
}

double summation_1d_sin_even (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	double res = 0 ;
	for (int i=1 ; i<nbr ; i++) 
		res += tab(i)*sin(2*i*xx) ;
	return res ;
}

double summation_1d_sin_odd (double xx, const Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	double res = 0 ;
	for (int i=0 ; i<nbr ; i++) 
		res += tab(i)*sin((2*i+1)*xx) ;
	return res ;
}


double summation_1d (int base, double xx, const Array<double>& tab) {

    static double (*summation_1d[NBR_MAX_BASE])(double, const Array<double>&) ;
    static bool premier_appel = true ;
    // Premier appel
    if (premier_appel) {
	premier_appel = false ;

	for (int i=0; i<NBR_MAX_BASE; i++)
	    summation_1d[i] = summation_1d_pasprevu ;

	summation_1d[CHEB] = summation_1d_cheb ;
	summation_1d[CHEB_EVEN] = summation_1d_cheb_even ;
	summation_1d[CHEB_ODD] = summation_1d_cheb_odd ;
	summation_1d[COSSIN] = summation_1d_cossin ;
	summation_1d[COS_EVEN] = summation_1d_cos_even ;
	summation_1d[COS_ODD] = summation_1d_cos_odd ;
	summation_1d[SIN_ODD] = summation_1d_sin_odd ;
	summation_1d[SIN_EVEN] = summation_1d_sin_even ;
	summation_1d[COS] = summation_1d_cos ;
	summation_1d[SIN] = summation_1d_sin ;
	summation_1d[LEG] = summation_1d_leg ;
	summation_1d[LEG_EVEN] = summation_1d_leg_even ;
	summation_1d[LEG_ODD] = summation_1d_leg_odd ;
	summation_1d[COSSIN_EVEN] = summation_1d_cossin_even ;	
	summation_1d[COSSIN_ODD] = summation_1d_cossin_odd ;
	}
	
        return summation_1d[base](xx, tab) ;
}

double Base_spectral::summation (const Point& num, const Array<double>& cf) const {

	Array<double>* courant = new Array<double>(cf) ;
	Dim_array nbr_coefs (cf.get_dimensions()) ;
	
	// Loop on dimensions (except the last):
	for (int d=0 ; d<ndim-1 ; d++) {
		int dim_output = ndim-1-d ;
		Dim_array nbr_output (dim_output) ;
		for (int k=0 ; k<dim_output ; k++)
			nbr_output.set(k) = nbr_coefs(k+d+1) ;
		Array<double> output (nbr_output) ;

		Index inout (nbr_output) ;
		Array<double> tab_1d (cf.get_size(d)) ;
		Index incourant (courant->get_dimensions()) ;
		
		// Loop on the points :
		bool loop = true ;
		while (loop) {
			int base = (*bases_1d[d])(inout) ;
			
			for (int k=0 ; k<dim_output ; k++)
				incourant.set(k+1) = inout(k) ;
			for (int k=0 ; k<cf.get_size(d) ; k++) {
				incourant.set(0) = k ;
				tab_1d.set(k) = (*courant)(incourant) ;
			}
			output.set(inout) = summation_1d (base, num(d+1), tab_1d) ;
			loop = inout.inc() ;
        }
		delete courant ;
		courant = new Array<double> (output) ;
	}
	
	// Last base :
	assert (courant->get_ndim()==1) ;
	assert (bases_1d[ndim-1]->get_ndim()==1) ;
	assert (bases_1d[ndim-1]->get_size(0)==1) ;
	double res = summation_1d ((*bases_1d[ndim-1])(0), num(ndim), *courant) ;
	delete courant ;
	
	return res ;
}}
