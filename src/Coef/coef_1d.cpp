/*
    Copyright 2017 Philippe Grandclement
              2018 Ludwig Jens Papenfort

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
#include "base_fftw.hpp"
#include "headcpp.hpp"
#include "array.cpp"
#include <fftw3.h>
#include "math.h"
#include <unordered_map>
#include <tuple>

namespace Kadath {
// fftw3 computes optimal algorithm once
// then we store it both in a map for later use
std::unordered_map<int, fftw_precomp_t> fftw_precomp_map;

// get or create buffer and plan
fftw_precomp_t& coef_1d_fftw(int n) {
  auto precomp_it = fftw_precomp_map.find(n);
  if(precomp_it == fftw_precomp_map.end())
    precomp_it = fftw_precomp_map.emplace(std::piecewise_construct,
                                          std::forward_as_tuple(n),
                                          std::forward_as_tuple(n,FFTW_R2HC)).first;
  return (*precomp_it).second;
}

void coef_i_1d (int, Array<double>&) ;

void coef_1d_pasprevu (Array<double>&) {
	cout << "Coef_1d not implemented." << endl ;
	abort() ;
}

void coef_1d_cheb (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	double fmoins0 = 0.5 * ( tab(0) - tab(nr-1));

  auto & fftw_data = coef_1d_fftw(nr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	for (int i=1 ; i<(nr-1)/2 ; i++) {
	     double fp = 0.5*(tab(i)+tab(nr-1-i)) ;
	     double fms = 0.5*(tab(i)-tab(nr-1-i))*sin(M_PI*i/(nr-1)) ;
	     tab_auxi[i] = fp + fms ;
	     tab_auxi[nr-1-i] = fp-fms ;
	     }
	tab_auxi[0] = 0.5*(tab(0) + tab(nr-1)) ;
	tab_auxi[(nr-1)/2] = tab((nr-1)/2) ;

	fftw_execute(p) ;

	// Coefficient pairs :
	tab.set(0) = tab_auxi[0] / (nr-1) ;
	for (int i=2; i<nr-1; i += 2)
	    tab.set(i) = 2*tab_auxi[(i/2)]/(nr-1) ;
	tab.set(nr-1) = tab_auxi[(nr-1)/2] / (nr-1) ;

	// Coefficients impairs :
    	tab.set(1) = 0 ;
    	double som = 0;
    	for (int  i = 3; i < nr; i += 2 ) {
	      tab.set(i) = tab(i-2) -4 * tab_auxi[nr-1-(i-1)/2]/(nr-1) ;
	      som += tab(i) ;
    	    }

// 2. Calcul de c_1 :
	    double c1 = - ( fmoins0 + som ) / ((nr-1)/2) ;

// 3. Coef. c_k avec k impair:
    	    tab.set(1) = c1 ;
    	    for (int i = 3; i < nr; i += 2 )
	    	tab.set(i) += c1 ;
}

void coef_1d_cheb_even (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	double fmoins0 = - 0.5 * ( tab(0) - tab(nr-1));

  auto & fftw_data = coef_1d_fftw(nr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	for (int i=1 ; i<(nr-1)/2 ; i++) {
	     double fp = 0.5*(tab(i)+tab(nr-1-i)) ;
	     double fms = 0.5*(-tab(i)+tab(nr-1-i))*sin(M_PI*i/(nr-1)) ;
	     tab_auxi[i] = fp + fms ;
	     tab_auxi[nr-1-i] = fp-fms ;
	     }
	tab_auxi[0] = 0.5*(tab(0) + tab(nr-1)) ;
	tab_auxi[(nr-1)/2] = tab((nr-1)/2) ;

	fftw_execute(p) ;

	// Coefficient pairs :
	tab.set(0) = tab_auxi[0] / (nr-1) ;
	for (int i=2; i<nr-1; i += 2)
	    tab.set(i) = 2*tab_auxi[(i/2)]/(nr-1) ;
	tab.set(nr-1) = tab_auxi[(nr-1)/2] / (nr-1) ;

	// Coefficients impairs :
    	tab.set(1) = 0 ;
    	double som = 0;
    	for (int  i = 3; i < nr; i += 2 ) {
	      tab.set(i) = tab(i-2) +4 * tab_auxi[nr-1-i/2]/(nr-1) ;
	      som += tab(i) ;
    	    }

// 2. Calcul de c_1 :
	    double c1 = ( fmoins0 - som ) / ((nr-1)/2) ;

// 3. Coef. c_k avec k impair:
    	    tab.set(1) = c1 ;
    	    for (int i = 3; i < nr; i += 2 )
	    	tab.set(i) += c1 ;
}

void coef_1d_cheb_odd (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nr = tab.get_size(0) ;

	double* cf = fftw_alloc_real(nr-1);
	for (int i=0 ; i<nr ; i++)
	    cf[i] = tab(i) * sin(M_PI/2.*i/(nr-1)) ;

	double fmoins0 = - 0.5 * ( cf[0] - cf[nr-1] );

  auto & fftw_data = coef_1d_fftw(nr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	for (int i=1 ; i<(nr-1)/2 ; i++) {
	     double fp = 0.5*(cf[i]+cf[nr-1-i]) ;
	     double fms = 0.5*(-cf[i]+cf[nr-1-i])*sin(M_PI*i/(nr-1)) ;
	     tab_auxi[i] = fp + fms ;
	     tab_auxi[nr-1-i] = fp-fms ;
	     }
	tab_auxi[0] = 0.5*(cf[0] + cf[nr-1]) ;
	tab_auxi[(nr-1)/2] = cf[(nr-1)/2] ;

	fftw_execute(p) ;

	// Coefficient pairs :
	cf[0] = tab_auxi[0] / (nr-1) ;
	for (int i=2; i<nr-1; i += 2)
	    cf[i] = 2*tab_auxi[(i/2)]/(nr-1) ;
	cf[nr-1] = tab_auxi[(nr-1)/2] / (nr-1) ;

	// Coefficients impairs :
    	cf[1] = 0 ;
    	double som = 0;
    	for (int  i = 3; i < nr; i += 2 ) {
	      cf[i] = cf[i-2] +4 * tab_auxi[nr-1-i/2]/(nr-1) ;
	      som += cf[i] ;
    	    }

// 2. Calcul de c_1 :
	    double c1 = ( fmoins0 - som ) / ((nr-1)/2) ;

// 3. Coef. c_k avec k impair:
    	    cf[1] = c1 ;
    	    for (int i = 3; i < nr; i += 2 )
	    	cf[i] += c1 ;

	cf[0] = 2* cf[0] ;
    	for (int i=1; i<nr-1; i++)
	    cf[i] = 2 * cf[i] - cf[i-1] ;
        cf[nr-1] = 0 ;

	for (int i=0 ; i<nr ; i++)
	    tab.set(i) = cf[i] ;

	delete [] cf ;
}

double coloc_leg(int,int) ;
double weight_leg(int,int) ;
double gamma_leg(int,int) ;
double leg(int,double) ;
double coloc_leg_parity(int,int) ;
double weight_leg_parity(int,int) ;
double gamma_leg_even(int,int) ;
double gamma_leg_odd(int,int) ;

void coef_1d_leg (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

	Array<double> res (nbr) ;
	res = 0 ;
	for (int i=0 ; i<nbr ; i++)
	    for (int j=0 ; j<nbr ; j++)
	     	res.set(i) += tab(j)*leg(i, coloc_leg(j,nbr))*weight_leg(j,nbr)/gamma_leg(i,nbr) ;
	tab = res ;
}

void coef_1d_leg_even (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

	Array<double> res (nbr) ;
	res = 0 ;
	for (int i=0 ; i<nbr ; i++)
	    for (int j=0 ; j<nbr ; j++)
	     	res.set(i) += tab(j)*leg(i*2, coloc_leg_parity(j,nbr))
				*weight_leg_parity(j,nbr)/gamma_leg_even(i,nbr) ;
	tab = res ;
}

void coef_1d_leg_odd (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

	Array<double> res (nbr) ;
	res = 0 ;
	for (int i=0 ; i<nbr-1 ; i++)
	    for (int j=0 ; j<nbr ; j++)
	     	res.set(i) += tab(j)*leg(i*2+1, coloc_leg_parity(j,nbr))
				*weight_leg_parity(j,nbr)/gamma_leg_odd(i,nbr) ;
	tab = res ;
}

void coef_1d_cossin (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	if (nbr>3) {

  auto & fftw_data = coef_1d_fftw(nbr-2);
	double* cf = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	for (int i=0 ; i<nbr-2 ; i++)
	     cf[i] = tab(i) ;
	fftw_execute(p) ;

	int index = 0 ;
	double* pcos = cf ;
	double* psin = cf+nbr-3 ;

	tab.set(index) = *pcos/double(nbr-2) ;
	index++ ; pcos ++ ;
	tab.set(index) = 0 ;
	index ++ ;

	for (int i=1 ; i<nbr/2 - 1; i++) {
	    tab.set(index) = 2.*(*pcos)/double(nbr-2) ; // the cosines
 	    index++ ; pcos ++ ;
 	    tab.set(index) = -2.*(*psin)/double(nbr-2) ; // the sines
  	    index ++ ; psin -- ;
 	}

	tab.set(index) = *pcos/double(nbr-2) ;
	index ++ ;
	tab.set(index) = 0 ;
	}
	else {
	  tab.set(1) = 0 ;
	  tab.set(2) = 0 ;
	}
}

void coef_1d_cos (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	// Symetrie taken into account
	double fmoins0 = 0.5 * ( tab(0) - tab(nbr-1) );

  auto & fftw_data = coef_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	for (int i = 1; i < (nbr-1)/2 ; i++ ) {
		double fp = 0.5 * ( tab(i) + tab(nbr-1-i) ) ;
		double fms = 0.5 * ( tab(i) - tab(nbr-1-i) ) * sin(M_PI*i/(nbr-1)) ;
		tab_auxi[i] = fp + fms ;
		tab_auxi[nbr-1-i] = fp - fms ;
    	    }

    	tab_auxi[0] = 0.5 * ( tab(0) + tab(nbr-1) );
    	tab_auxi[(nbr-1)/2] = tab((nbr-1)/2);

	fftw_execute(p) ;

	tab.set(0) = tab_auxi[0] / (nbr-1) ;
    	for (int i=2; i<nbr-1; i += 2 )
		tab.set(i) = 2*tab_auxi[i/2]/(nbr-1) ;
	tab.set(nbr-1) = tab_auxi[(nbr-1)/2] / (nbr-1) ;

	tab.set(1) = 0 ;
    	double som = 0;
    	for (int i = 3; i < nbr; i += 2 ) {
		tab.set(i) = tab(i-2) + 4 * tab_auxi[nbr-1 - i/2]/(nbr-1) ;
	    	som += tab(i) ;
    	}

// 2. Calcul de c_1 :
    	    double c1 = ( fmoins0 - som ) / ((nbr-1)/2) ;

// 3. Coef. c_k avec k impair:
    	    tab.set(1) = c1 ;
    	    for (int i = 3; i < nbr; i += 2 )
	    	tab.set(i) += c1 ;
}

void coef_1d_sin (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	// Symetrie taken into account
  auto & fftw_data = coef_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	for (int i = 1; i < (nbr-1)/2 ; i++ ) {
		double fp = 0.5 * ( tab(i) + tab(nbr-1-i) ) * sin(M_PI*i/(nbr-1)) ;
		double fms = 0.5 * ( tab(i) - tab(nbr-1-i) )  ;
		tab_auxi[i] = fp + fms ;
		tab_auxi[nbr-1-i] = fp - fms ;
    	    }

    	tab_auxi[0] = 0.5 * ( tab(0) + tab(nbr-1) );
    	tab_auxi[(nbr-1)/2] = tab((nbr-1)/2);

	fftw_execute(p) ;

	tab.set(0) = 0 ;
    	for (int i=2; i<nbr-1; i += 2 )
		tab.set(i) = -2*tab_auxi[nbr-1-i/2]/(nbr-1) ;
	tab.set(nbr-1) = 0 ;

	tab.set(1) = 2*tab_auxi[0]/(nbr-1) ;
    	for (int i = 3; i < nbr; i += 2 )
		tab.set(i) = tab(i-2) + 4 * tab_auxi[i/2]/(nbr-1) ;
}

void coef_1d_cos_even (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	if (nbr>3) {
	// Symetrie taken into account
	double fmoins0 = 0.5 * ( tab(0) - tab(nbr-1) );

  auto & fftw_data = coef_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	for (int i = 1; i < (nbr-1)/2 ; i++ ) {
		double fp = 0.5 * ( tab(i) + tab(nbr-1-i) ) ;
		double fms = 0.5 * ( tab(i) - tab(nbr-1-i) ) * sin(M_PI*i/(nbr-1)) ;
		tab_auxi[i] = fp + fms ;
		tab_auxi[nbr-1-i] = fp - fms ;
    	    }

    	tab_auxi[0] = 0.5 * ( tab(0) + tab(nbr-1) );
    	tab_auxi[(nbr-1)/2] = tab((nbr-1)/2);

	fftw_execute(p) ;

	tab.set(0) = tab_auxi[0] / (nbr-1) ;
    	for (int i=2; i<nbr-1; i += 2 )
		tab.set(i) = 2*tab_auxi[i/2]/(nbr-1) ;
	tab.set(nbr-1) = tab_auxi[(nbr-1)/2] / (nbr-1) ;

	tab.set(1) = 0 ;
    	double som = 0;
    	for (int i = 3; i < nbr; i += 2 ) {
		tab.set(i) = tab(i-2) + 4 * tab_auxi[nbr-1 - i/2]/(nbr-1) ;
	    	som += tab(i) ;
    	}

// 2. Calcul de c_1 :
    	    double c1 = ( fmoins0 - som ) / ((nbr-1)/2) ;

// 3. Coef. c_k avec k impair:
    	    tab.set(1) = c1 ;
    	    for (int i = 3; i < nbr; i += 2 )
	    	tab.set(i) += c1 ;
	}
}

void coef_1d_cos_odd (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	if (nbr>3) {
	// Symetrie taken into account
	double* cf = fftw_alloc_real(nbr);
	for (int i=0 ; i<nbr-1 ; i++)
	    cf[i] = tab(i)*sin(M_PI*(nbr-1-i)/2./(nbr-1)) ;
	cf[nbr-1] = 0 ;
	double fmoins0 = 0.5 * (cf[0] - cf[nbr-1] );

  auto & fftw_data = coef_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	  for (int i = 1; i < (nbr-1)/2 ; i++ ) {
		double fp = 0.5 * ( cf[i] + cf[nbr-1-i] ) ;
		double fms = 0.5 * ( cf[i] - cf[nbr-1-i] ) * sin(M_PI*i/(nbr-1)) ;
		tab_auxi[i] = fp + fms ;
		tab_auxi[nbr-1-i] = fp - fms ;
    	    }

    	tab_auxi[0] = 0.5 * ( cf[0] + cf[nbr-1] );
    	tab_auxi[(nbr-1)/2] = cf[(nbr-1)/2];

	fftw_execute(p) ;

	cf[0] = tab_auxi[0] / (nbr-1) ;
    	for (int i=2; i<nbr-1; i += 2 )
		cf[i] = 2*tab_auxi[i/2]/(nbr-1) ;
	cf[nbr-1] = tab_auxi[(nbr-1)/2] ;

	cf[1] = 0 ;
    	double som = 0;
    	for (int i = 3; i < nbr; i += 2 ) {
		cf[i] = cf[i-2] + 4 * tab_auxi[nbr-1 - i/2]/(nbr-1) ;
	    	som += cf[i] ;
    	}

// 2. Calcul de c_1 :
    	    double c1 = ( fmoins0 - som ) / ((nbr-1)/2) ;

// 3. Coef. c_k avec k impair:
    	    cf[1] = c1 ;
    	    for (int i = 3; i < nbr; i += 2 )
	    	cf[i] += c1 ;

	cf[0] = 2* cf[0] ;
    	for (int i=1; i<nbr-1; i++)
		cf[i] = 2 * cf[i] - cf[i-1] ;
    	cf[nbr-1] = 0 ;
	for (int i=0 ; i<nbr ; i++)
		tab.set(i) = cf[i] ;

	delete [] cf ;
	}
}

void coef_1d_sin_even (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	if (nbr>3) {
	// Symetrie taken into account
  auto & fftw_data = coef_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	for (int i = 1; i < (nbr-1)/2 ; i++ ) {
		double fp = 0.5 * ( tab(i) + tab(nbr-1-i) ) * sin(M_PI*i/(nbr-1));
		double fms = 0.5 * ( tab(i) - tab(nbr-1-i) )  ;
		tab_auxi[i] = fp + fms ;
		tab_auxi[nbr-1-i] = fp - fms ;
    	    }

    	tab_auxi[0] = 0.5 * ( tab(0) -tab(nbr-1) );
    	tab_auxi[(nbr-1)/2] = tab((nbr-1)/2);

	fftw_execute(p) ;

	tab.set(0) = 0. ;
    	for (int i=2; i<nbr-1; i += 2 )
		tab.set(i) = -2*tab_auxi[nbr-1-i/2]/(nbr-1) ;
	tab.set(nbr-1) = 0 ;

	tab.set(1) = 2*tab_auxi[0]/(nbr-1) ;
    	for (int i = 3; i < nbr; i += 2 )
		tab.set(i) = tab(i-2) + 4 * tab_auxi[i/2]/(nbr-1) ;
	}
}

void coef_1d_sin_odd (Array<double>& tab) {
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	if (nbr>3) {
	// Symetrie taken into account
	double* cf = new double[nbr] ;
	cf[0] = 0 ;
	for (int i=1 ; i<nbr ; i++)
	    cf[i] = tab(i)*sin(M_PI*i/2./(nbr-1)) ;
	double fmoins0 = 0.5 * (cf[0] - cf[nbr-1] );

  auto & fftw_data = coef_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	  for (int i = 1; i < (nbr-1)/2 ; i++ ) {
		double fp = 0.5 * ( cf[i] + cf[nbr-1-i] ) ;
		double fms = 0.5 * ( cf[i] - cf[nbr-1-i] ) * sin(M_PI*i/(nbr-1)) ;
		tab_auxi[i] = fp + fms ;
		tab_auxi[nbr-1-i] = fp - fms ;
    	    }

    	tab_auxi[0] = 0.5 * ( cf[0] + cf[nbr-1] );
    	tab_auxi[(nbr-1)/2] = cf[(nbr-1)/2];

	fftw_execute(p) ;

	cf[0] = tab_auxi[0] / (nbr-1) ;
    	for (int i=2; i<nbr-1; i += 2 )
		cf[i] = 2*tab_auxi[i/2]/(nbr-1) ;
	cf[nbr-1] = tab_auxi[(nbr-1)/2] / (nbr-1) ;

	cf[1] = 0 ;
    	double som = 0;
    	for (int i = 3; i < nbr; i += 2 ) {
		cf[i] = cf[i-2] + 4 * tab_auxi[nbr-1 - i/2]/(nbr-1) ;
	    	som += cf[i] ;
    	}

// 2. Calcul de c_1 :
    	    double c1 = ( fmoins0 - som ) / ((nbr-1)/2) ;

// 3. Coef. c_k avec k impair:
    	    cf[1] = c1 ;
    	    for (int i = 3; i < nbr; i += 2 )
	    	cf[i] += c1 ;

	cf[0] = 2* cf[0] ;
    	for (int i=1; i<nbr-1; i++)
		cf[i] = 2 * cf[i] + cf[i-1] ;
    	cf[nbr-1] = 0 ;
	for (int i=0 ; i<nbr ; i++)
		tab.set(i) = cf[i] ;

	delete [] cf ;
	}
}


void coef_1d_cossin_even (Array<double>& tab) {

	// Copy values in double-size array :
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	Array<double> tab2(2*nbr-2) ;
	for (int i=0 ; i<nbr-2 ; i++) {
	  tab2.set(i) = tab(i) ;
	  tab2.set(i+nbr-2) = tab(i) ;
	}

	coef_1d_cossin (tab2) ;

	int conte = 0 ;
	for (int k=0 ; k<nbr-1 ; k+=2) {
	  tab.set(k) = tab2(conte) ;
	  tab.set(k+1) = tab2(conte+1) ;
	  conte += 4 ;
	}
}


void coef_1d_cossin_odd (Array<double>& tab) {
	// Copy values in double-size array :
	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	Array<double> tab2(2*nbr-2) ;
	tab2 = 0 ;
	for (int i=0 ; i<nbr-2 ; i++) {
	  tab2.set(i) = tab(i) ;
	  tab2.set(i+nbr-2) = -tab(i) ;
	}
	coef_1d_cossin (tab2) ;

	int conte = 2 ;
	for (int k=0 ; k<nbr-3 ; k+=2) {
	  tab.set(k) = tab2(conte) ;
	  tab.set(k+1) = tab2(conte+1) ;
	  conte += 4 ;
	}
	tab.set(nbr-2) = 0. ;
	tab.set(nbr-1) = 0. ;
}

void coef_1d (int base, Array<double>& tab) {
    static void (*coef_1d[NBR_MAX_BASE])(Array<double>&) ;
    static bool premier_appel = true ;

    // Premier appel
    if (premier_appel) {
	premier_appel = false ;

	for (int i=0; i<NBR_MAX_BASE; i++)
	    coef_1d[i] = coef_1d_pasprevu ;

	coef_1d[CHEB] = coef_1d_cheb ;
	coef_1d[CHEB_EVEN] = coef_1d_cheb_even ;
	coef_1d[CHEB_ODD] = coef_1d_cheb_odd ;
	coef_1d[COSSIN] = coef_1d_cossin ;
	coef_1d[COS_EVEN] = coef_1d_cos_even ;
	coef_1d[COS_ODD] = coef_1d_cos_odd ;
	coef_1d[SIN_EVEN] = coef_1d_sin_even ;
	coef_1d[SIN_ODD] = coef_1d_sin_odd ;
	coef_1d[COS] = coef_1d_cos ;
	coef_1d[SIN] = coef_1d_sin ;
	coef_1d[LEG] = coef_1d_leg ;
	coef_1d[LEG_EVEN] = coef_1d_leg_even ;
	coef_1d[LEG_ODD] = coef_1d_leg_odd;
	coef_1d[COSSIN_EVEN] = coef_1d_cossin_even ;
	coef_1d[COSSIN_ODD] = coef_1d_cossin_odd ;
	}

        coef_1d[base](tab) ;
}
}

