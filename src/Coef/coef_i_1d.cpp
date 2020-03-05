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
#include <unordered_map>
#include <tuple>

namespace Kadath {
// fftw3 computes optimal algorithm once
// then we store it both in a map for later use
std::unordered_map<int, fftw_precomp_t> fftw_precomp_map_i;

// get or create buffer and plan
fftw_precomp_t& coef_i_1d_fftw(int n) {
  auto precomp_it = fftw_precomp_map_i.find(n);
  if(precomp_it == fftw_precomp_map_i.end())
    precomp_it = fftw_precomp_map_i.emplace(std::piecewise_construct,
                                            std::forward_as_tuple(n),
                                            std::forward_as_tuple(n,FFTW_HC2R)).first;
  return (*precomp_it).second;
}

void coef_i_1d_pasprevu (Array<double>&) {
	cout << "Coef_1d not implemented." << endl ;
	abort() ;
}

void coef_i_1d_cheb (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

	if (nbr>3) {

  auto & fftw_data = coef_i_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	double* cf = new double[nbr];

	double c1 = tab(1) ;
	double somme = 0 ;
	cf[1] = 0 ;
	for (int i=3 ; i<nbr ; i+=2) {
		cf[i] = tab(i) - c1 ;
		somme += cf[i] ;
	}
	double fmoins0 = -(nbr-1)/2 *c1- somme ;
	for (int i=3 ; i<nbr ; i+=2)
		tab_auxi[nbr-1-i/2] = -0.25*(cf[i]-cf[i-2]) ;
	tab_auxi[0] = tab(0) ;
	for (int i=1 ; i<(nbr-1)/2 ; i++)
		tab_auxi[i] = 0.5*tab(2*i) ;
	tab_auxi[(nbr-1)/2] = tab(nbr-1) ;

	fftw_data.execute() ;

	for (int i=1 ; i<(nbr-1)/2 ; i++) {
		double fp = 0.5*(tab_auxi[i]+tab_auxi[nbr-1-i]) ;
		double fm = 0.5*(tab_auxi[i]-tab_auxi[nbr-1-i])/sin(M_PI*i/(nbr-1)) ;
		tab.set(i) = fp+fm ;
		tab.set(nbr-i-1) = fp -fm ;
	}
	tab.set(0) = tab_auxi[0] + fmoins0 ;
	tab.set(nbr-1) = tab_auxi[0] - fmoins0 ;
	tab.set((nbr-1)/2) = tab_auxi[(nbr-1)/2] ;

	delete [] cf ;
	}
}

void coef_i_1d_cheb_even (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

  auto & fftw_data = coef_i_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	double* cf = new double[nbr] ;

	double c1 = tab(1) ;
	double somme = 0 ;
	cf[1] = 0 ;
	for (int i=3 ; i<nbr ; i+=2) {
		cf[i] = tab(i) - c1 ;
		somme += cf[i] ;
	}
	double fmoins0 = (nbr-1)/2 *c1+somme ;
	for (int i=3 ; i<nbr ; i+=2)
		tab_auxi[nbr-1-i/2] = 0.25*(cf[i]-cf[i-2]) ;
	tab_auxi[0] = tab(0) ;
	for (int i=1 ; i<(nbr-1)/2 ; i++)
		tab_auxi[i] = 0.5*tab(2*i) ;
	tab_auxi[(nbr-1)/2] = tab(nbr-1) ;

	fftw_data.execute() ;

	for (int i=1 ; i<(nbr-1)/2 ; i++) {
		double fp = 0.5*(tab_auxi[i]+tab_auxi[nbr-1-i]) ;
		double fm = 0.5*(tab_auxi[i]-tab_auxi[nbr-1-i])/sin(M_PI*i/(nbr-1)) ;
		tab.set(nbr-1-i) = fp+fm ;
		tab.set(i) = fp -fm ;
	}
	tab.set(0) = tab_auxi[0] - fmoins0 ;
	tab.set(nbr-1) = tab_auxi[0] + fmoins0 ;
	tab.set((nbr-1)/2) = tab_auxi[(nbr-1)/2] ;

	delete [] cf ;
}

void coef_i_1d_cheb_odd (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

  auto & fftw_data = coef_i_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	double* ti = new double [nbr] ;
	double* cf = new double [nbr] ;

	ti[0] = 0.5*tab(0) ;
	for (int i=1 ; i<nbr-1 ; i++)
		ti[i] = 0.5*(tab(i)+tab(i-1)) ;
	ti[nbr-1] = 0.5*tab(nbr-2) ;

	double c1 = ti[1] ;
	double somme = 0 ;
	cf[1] = 0 ;
	for (int i=3 ; i<nbr ; i+=2) {
		cf[i] = ti[i] - c1 ;
		somme += cf[i] ;
	}
	double fmoins0 = (nbr-1)/2 *c1+ somme ;
	for (int i=3 ; i<nbr ; i+=2)
		tab_auxi[nbr-1-i/2] = 0.25*(cf[i]-cf[i-2]) ;
	tab_auxi[0] = ti[0] ;
	for (int i=1 ; i<(nbr-1)/2 ; i++)
		tab_auxi[i] = 0.5*ti[2*i] ;
	tab_auxi[(nbr-1)/2]=ti[nbr-1] ;

	fftw_data.execute() ;

	for (int i=1 ; i<(nbr-1)/2 ; i++) {
		double fp = 0.5*(tab_auxi[i]+tab_auxi[nbr-1-i]) ;
		double fm = 0.5*(tab_auxi[i]-tab_auxi[nbr-1-i])/sin(M_PI*i/(nbr-1)) ;
		tab.set(nbr-1-i) = (fp+fm)/sin(M_PI*(nbr-1-i)/2/(nbr-1)) ;
		tab.set(i) = (fp -fm)/sin(M_PI*i/2/(nbr-1)) ;
	}
	tab.set(0) = 0 ;
	tab.set(nbr-1) = tab_auxi[0] + fmoins0 ;
	tab.set((nbr-1)/2) = tab_auxi[(nbr-1)/2]/sin(M_PI*(nbr-1)/4/(nbr-1)) ;

	delete [] ti ;
	delete [] cf ;
}

double coloc_leg(int,int) ;
double summation_1d_leg(double xx, const Array<double>&) ;
double coloc_leg_parity(int,int) ;
double summation_1d_leg_even(double xx, const Array<double>&) ;
double summation_1d_leg_odd(double xx, const Array<double>&) ;

void coef_i_1d_leg (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

	Array<double> res(nbr) ;
	for (int i=0 ; i<nbr ; i++)
		res.set(i) = summation_1d_leg(coloc_leg(i,nbr), tab) ;
	tab=res ;
}

void coef_i_1d_leg_even (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

	Array<double> res(nbr) ;
	for (int i=0 ; i<nbr ; i++)
		res.set(i) = summation_1d_leg_even(coloc_leg_parity(i,nbr), tab) ;
	tab=res ;
}

void coef_i_1d_leg_odd (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

	Array<double> res(nbr) ;
	for (int i=0 ; i<nbr ; i++)
		res.set(i) = summation_1d_leg_odd(coloc_leg_parity(i,nbr), tab) ;
	tab=res ;
}

void coef_i_1d_cossin (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	int np = nbr-2 ;
	if (np>1) {

  auto & fftw_data = coef_i_1d_fftw(np);
	double* cf = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	cf[0] = tab(0) ;
	for (int i=1 ; i<np/2 ; i++) {
	     cf[i] = 0.5*tab(2*i) ;
	     cf[np-i] = -0.5*tab(2*i+1) ;
	     }
	cf[np/2] = tab(np) ;
	fftw_data.execute() ;

	for (int i=0 ; i<np ; i++)
		tab.set(i) = cf[i] ;
	}
}

void coef_i_1d_cos (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

  auto & fftw_data = coef_i_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	double* cf = new double [nbr] ;

	double c1 = tab(1) ;
	double somme = 0 ;
	cf[1] = 0 ;
	for (int i=3 ; i<nbr ; i+=2) {
		cf[i] = tab(i) - c1 ;
		somme += cf[i] ;
	}
	double fmoins0 = (nbr-1)/2 *c1+ somme ;
	for (int i=3 ; i<nbr ; i+=2)
		tab_auxi[nbr-1-i/2] = 0.25*(cf[i]-cf[i-2]) ;
	tab_auxi[0] = tab(0) ;
	for (int i=1 ; i<(nbr-1)/2 ; i++)
		tab_auxi[i] = 0.5*tab(2*i) ;
	tab_auxi[(nbr-1)/2] = tab(nbr-1) ;

	fftw_data.execute() ;

	for (int i=1 ; i<(nbr-1)/2 ; i++) {
		double fp = 0.5*(tab_auxi[i]+tab_auxi[nbr-1-i]) ;
		double fm = 0.5*(tab_auxi[i]-tab_auxi[nbr-1-i])/sin(M_PI*i/(nbr-1)) ;
		tab.set(i) = fp+fm ;
		tab.set(nbr-i-1) = fp -fm ;
	}
	tab.set(0) = tab_auxi[0] + fmoins0 ;
	tab.set(nbr-1) = tab_auxi[0] - fmoins0 ;
	tab.set((nbr-1)/2) = tab_auxi[(nbr-1)/2] ;

	delete [] cf;
}

void coef_i_1d_sin (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;

  auto & fftw_data = coef_i_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	for (int i=2 ; i<nbr-1 ; i+=2)
	    tab_auxi[nbr-1-i/2] = -0.5*tab(i) ;
	tab_auxi[0] = 0.5*tab(1) ;
	for (int i=3 ; i<nbr ; i+=2)
		tab_auxi[i/2] = 0.25*(tab(i)-tab(i-2)) ;
	tab_auxi[(nbr-1)/2] = -0.5*tab(nbr-2) ;

	fftw_data.execute() ;

	for (int i=1 ; i<(nbr-1)/2 ; i++) {
		double fp = 0.5*(tab_auxi[i]+tab_auxi[nbr-1-i])/sin(M_PI*i/(nbr-1)) ;
		double fm = 0.5*(tab_auxi[i]-tab_auxi[nbr-1-i]) ;
		tab.set(i) = fp+fm ;
		tab.set(nbr-i-1) = fp -fm ;
	}
	tab.set(0) = 0 ;
	tab.set(nbr-1) =-2* tab_auxi[0] ;
	tab.set((nbr-1)/2) = tab_auxi[(nbr-1)/2] ;
}

void coef_i_1d_cos_even (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	if (nbr>3) {
	double* cf = new double [nbr] ;

  auto & fftw_data = coef_i_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	double c1 = tab(1) ;
	double somme = 0 ;
	cf[1] = 0 ;
	for (int i=3 ; i<nbr ; i+=2) {
		cf[i] = tab(i) - c1 ;
		somme += cf[i] ;
	}
	double fmoins0 = (nbr-1)/2 *c1+ somme ;
	for (int i=3 ; i<nbr ; i+=2)
		tab_auxi[nbr-1-i/2] = 0.25*(cf[i]-cf[i-2]) ;
	tab_auxi[0] = tab(0) ;
	for (int i=1 ; i<(nbr-1)/2 ; i++)
		tab_auxi[i] = 0.5*tab(2*i) ;
	tab_auxi[(nbr-1)/2] = tab(nbr-1) ;

	fftw_data.execute() ;

	for (int i=1 ; i<(nbr-1)/2 ; i++) {
		double fp = 0.5*(tab_auxi[i]+tab_auxi[nbr-1-i]) ;
		double fm = 0.5*(tab_auxi[i]-tab_auxi[nbr-1-i])/sin(M_PI*i/(nbr-1)) ;
		tab.set(i) = fp+fm ;
		tab.set(nbr-i-1) = fp -fm ;
	}
	tab.set(0) = tab_auxi[0] + fmoins0 ;
	tab.set(nbr-1) = tab_auxi[0] - fmoins0 ;
	tab.set((nbr-1)/2) = tab_auxi[(nbr-1)/2] ;

	delete [] cf;
	}
}

void coef_i_1d_cos_odd (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	if (nbr>3) {
	double* ti = new double[nbr] ;
	double* cf = new double[nbr] ;

  auto & fftw_data = coef_i_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	ti[0] = 0.5*tab(0) ;
	for (int i=1 ; i<nbr-1 ; i++)
		ti[i] = 0.5*(tab(i)+tab(i-1)) ;
	ti[nbr-1] = 0.5*tab(nbr-2) ;

	double c1 = ti[1] ;
	double somme = 0 ;
	cf[1] = 0 ;
	for (int i=3 ; i<nbr ; i+=2) {
		cf[i] = ti[i] - c1 ;
		somme += cf[i] ;
	}
	double fmoins0 = (nbr-1)/2 *c1+ somme ;
	for (int i=3 ; i<nbr ; i+=2)
		tab_auxi[nbr-1-i/2] = 0.25*(cf[i]-cf[i-2]) ;
	tab_auxi[0] = ti[0] ;
	for (int i=1 ; i<(nbr-1)/2 ; i++)
		tab_auxi[i] = 0.5*ti[2*i] ;
	tab_auxi[(nbr-1)/2] = ti[nbr-1] ;

	fftw_data.execute() ;

	for (int i=1 ; i<(nbr-1)/2 ; i++) {
		double fp = 0.5*(tab_auxi[i]+tab_auxi[nbr-1-i]) ;
		double fm = 0.5*(tab_auxi[i]-tab_auxi[nbr-1-i])/sin(M_PI*i/(nbr-1)) ;
		tab.set(i) = (fp+fm)/sin(M_PI*(nbr-1-i)/2/(nbr-1)) ;
		tab.set(nbr-i-1) = (fp -fm)/sin(M_PI*i/2/(nbr-1)) ;
	}
	tab.set(0) = tab_auxi[0] + fmoins0 ;
	tab.set(nbr-1) = 0 ;
	tab.set((nbr-1)/2) = tab_auxi[(nbr-1)/2]/sin(M_PI*(nbr-1)/4/(nbr-1)) ;

	delete [] ti ;
	delete [] cf ;
	}
}

void coef_i_1d_sin_even (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	if (nbr>3) {

  auto & fftw_data = coef_i_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	for (int i=2 ; i<nbr-1 ; i+=2)
	    tab_auxi[nbr-1-i/2] = -0.5*tab(i) ;
	tab_auxi[0] = 0.5*tab(1) ;
	for (int i=1 ; i<(nbr-1)/2 ; i++)
		tab_auxi[i] = 0.25*(tab(2*i+1)-tab(2*i-1)) ;
	tab_auxi[(nbr-1)/2] = -0.5*tab(nbr-2) ;

	fftw_data.execute() ;

	for (int i=1 ; i<(nbr-1)/2 ; i++) {
		double fp = 0.5*(tab_auxi[i]+tab_auxi[nbr-1-i])/sin(M_PI*i/(nbr-1)) ;
		double fm = 0.5*(tab_auxi[i]-tab_auxi[nbr-1-i]) ;
		tab.set(i) = fp+fm ;
		tab.set(nbr-i-1) = fp -fm ;
	}
	tab.set(0) = 0 ;
	tab.set(nbr-1) =-2* tab_auxi[0] ;
	tab.set((nbr-1)/2) = tab_auxi[(nbr-1)/2] ;
	}
}

void coef_i_1d_sin_odd (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	if (nbr>3) {
	double* ti = new double[nbr] ;
	double* cf = new double[nbr] ;

  auto & fftw_data = coef_i_1d_fftw(nbr-1);
	double* tab_auxi = fftw_data.buffer;
	fftw_plan p = fftw_data.plan;

	ti[0] = 0.5*tab(0) ;
	for (int i=1 ; i<nbr-1 ; i++)
		ti[i] = 0.5*(tab(i)-tab(i-1)) ;
	ti[nbr-1] = -0.5*tab(nbr-2) ;

	double c1 = ti[1] ;
	double somme = 0 ;
	cf[1] = 0 ;
	for (int i=3 ; i<nbr ; i+=2) {
		cf[i] = ti[i] - c1 ;
		somme += cf[i] ;
	}
	double fmoins0 = (nbr-1)/2 *c1+ somme ;
	for (int i=3 ; i<nbr ; i+=2)
		tab_auxi[nbr-1-i/2] = 0.25*(cf[i]-cf[i-2]) ;
	tab_auxi[0] = ti[0] ;
	for (int i=1 ; i<(nbr-1)/2 ; i++)
		tab_auxi[i] = 0.5*ti[2*i] ;
	tab_auxi[(nbr-1)/2] = ti[nbr-1] ;

	fftw_data.execute() ;

	for (int i=1 ; i<(nbr-1)/2 ; i++) {
		double fp = 0.5*(tab_auxi[i]+tab_auxi[nbr-1-i]) ;
		double fm = 0.5*(tab_auxi[i]-tab_auxi[nbr-1-i])/sin(M_PI*i/(nbr-1)) ;
		tab.set(i) = (fp+fm)/sin(M_PI*i/2/(nbr-1)) ;
		tab.set(nbr-i-1) = (fp -fm)/sin(M_PI*(nbr-1-i)/2/(nbr-1)) ;
	}
	tab.set(0) = 0 ;
	tab.set(nbr-1) = tab_auxi[0] - fmoins0 ;
	tab.set((nbr-1)/2) = tab_auxi[(nbr-1)/2]/sin(M_PI*(nbr-1)/4/(nbr-1)) ;

	delete [] ti;
	delete [] cf ;
	}
}


void coef_i_1d_cossin_even (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	// Double-sized array :
	Array<double> tab2 (nbr*2-2) ;
	int conte = 0 ;
	for (int i=0 ; i<nbr-2 ; i+=2) {
	    tab2.set(conte) = tab(i) ;
	    tab2.set(conte+1) = tab(i+1) ;
	    tab2.set(conte+2) = 0. ;
	    tab2.set(conte+3) = 0. ;
	    conte += 4 ;
	}
	tab2.set(conte) = tab(nbr-2) ;
	tab2.set(conte+1) = tab(nbr-1) ;

	coef_i_1d_cossin(tab2) ;
	for (int i=0 ; i<nbr-2 ; i++)
	  tab.set(i) = tab2(i) ;
	tab.set(nbr-2) = 0 ;
	tab.set(nbr-1) = 0 ;
}


void coef_i_1d_cossin_odd (Array<double>& tab) {

	assert (tab.get_ndim()==1) ;
	int nbr = tab.get_size(0) ;
	// Double-sized array :
	Array<double> tab2 (nbr*2-2) ;
	int conte = 0 ;
	for (int i=0 ; i<nbr-2 ; i+=2) {
	    tab2.set(conte) = 0. ;
	    tab2.set(conte+1) = 0. ;
	    tab2.set(conte+2) = tab(i) ;
	    tab2.set(conte+3) = tab(i+1) ;
	    conte += 4 ;
	}
	tab2.set(conte) = tab(nbr-2) ;
	tab2.set(conte+1) = tab(nbr-1) ;

	coef_i_1d_cossin(tab2) ;
	for (int i=0 ; i<nbr-2 ; i++)
	  tab.set(i) = tab2(i) ;
      	tab.set(nbr-2) = 0 ;
	tab.set(nbr-1) = 0 ;
}


void coef_i_1d (int base, Array<double>& tab) {
    static void (*coef_i_1d[NBR_MAX_BASE])(Array<double>&) ;
    static bool premier_appel = true ;

    // Premier appel
    if (premier_appel) {
	premier_appel = false ;

	for (int i=0; i<NBR_MAX_BASE; i++)
	    coef_i_1d[i] = coef_i_1d_pasprevu ;

	coef_i_1d[CHEB] = coef_i_1d_cheb ;
	coef_i_1d[CHEB_EVEN] = coef_i_1d_cheb_even ;
	coef_i_1d[CHEB_ODD] = coef_i_1d_cheb_odd ;
	coef_i_1d[COSSIN] = coef_i_1d_cossin ;
	coef_i_1d[COS_EVEN] = coef_i_1d_cos_even ;
	coef_i_1d[COS_ODD] = coef_i_1d_cos_odd ;
	coef_i_1d[SIN_EVEN] = coef_i_1d_sin_even ;
	coef_i_1d[SIN_ODD] = coef_i_1d_sin_odd ;
	coef_i_1d[COS] = coef_i_1d_cos ;
	coef_i_1d[SIN] = coef_i_1d_sin ;
	coef_i_1d[LEG] = coef_i_1d_leg ;
	coef_i_1d[LEG_EVEN] = coef_i_1d_leg_even ;
	coef_i_1d[LEG_ODD] = coef_i_1d_leg_odd;
	coef_i_1d[COSSIN_EVEN] = coef_i_1d_cossin_even ;
	coef_i_1d[COSSIN_ODD] = coef_i_1d_cossin_odd ;
	}

        coef_i_1d[base](tab) ;
}
}

