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

#include "bispheric.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "param.hpp"
#include "utilities.hpp"
namespace Kadath {
double eta_lim_chi (double chi, double rext, double a, double eta_c) ;
double chi_lim_eta (double chi, double rext, double a, double chi_c) ;

double zerosec(double (*f)(double, const Param&), const Param& parf, 
    double x1, double x2, double precis, int nitermax, int& niter) ;
    

double func_a (double aa, const Param& par) {
	double r1 = par.get_double(0) ;
	double r2 = par.get_double(1) ;
	double d = par.get_double(2) ;
	return (sqrt(aa*aa+r1*r1)+sqrt(aa*aa+r2*r2)-d) ;
}

Space_bispheric::Space_bispheric (int ttype, double distance, double r1, double r2, double rext, int nr) {

    // Verif :
    ndim = 3 ;
    
    nbr_domains = 8 ;
    
    ndom_minus = 1 ;
    ndom_plus = 1 ;
    nshells = 0 ;
    
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
  
    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;
    
    // Bispheric :
    // Computation of aa
    Param par_a ;
    par_a.add_double(r1,0) ;
    par_a.add_double(r2,1) ;
    par_a.add_double(distance,2) ;
    double a_min = 0 ;
    double a_max = distance/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_a, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/r2) ;
    double eta_minus = -asinh(aa/r1) ;
    
    double chi_c = 2*atan(aa/rext) ;
    double eta_c = log((1+rext/aa)/(rext/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;
   
    // The spheres  
    Point center (ndim) ;
    center.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ; center.set(2) = 0 ; center.set(3) = 0 ;
    a_minus = aa*cosh(eta_minus)/sinh(eta_minus) ;
    domains[0] = new Domain_nucleus(0, ttype, r1, center, res) ;
    center.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    a_plus = aa*cosh(eta_plus)/sinh(eta_plus) ;
    domains[1] = new Domain_nucleus(1, ttype, r2, center, res) ;

    // Bispherical domains antitrigo order: 
    domains[2] = new Domain_bispheric_chi_first(2, ttype, aa, eta_minus, rext, chi_lim, res_bi) ;
    domains[3] = new Domain_bispheric_rect(3, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[4] = new Domain_bispheric_eta_first(4, ttype, aa, rext, -eta_lim, eta_lim, res_bi) ;
    domains[5] = new Domain_bispheric_rect(5, ttype, aa, rext, eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[6] = new Domain_bispheric_chi_first(6, ttype, aa, eta_plus, rext, chi_lim, res_bi) ;
    
    // ZEC
    center.set(1) = 0 ;
    domains[7] = new Domain_compact (7, ttype,rext, center, res) ;;
}

Space_bispheric::Space_bispheric (int ttype, double distance, double r1, double r2, int nn, const Array<double>& rr, int nr) {

    // Verif :
    ndim = 3 ;
    
     
    ndom_minus = 1 ;
    ndom_plus = 1 ;
    nshells = nn ;
    
    nbr_domains = 8 + nshells ;
  
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
  
    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;
    
    // Bispheric :
    // Computation of aa
    Param par_a ;
    par_a.add_double(r1,0) ;
    par_a.add_double(r2,1) ;
    par_a.add_double(distance,2) ;
    double a_min = 0 ;
    double a_max = distance/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_a, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/r2) ;
    double eta_minus = -asinh(aa/r1) ;

    double chi_c = 2*atan(aa/rr(0)) ;
    double eta_c = log((1+rr(0)/aa)/(rr(0)/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rr(0), aa, chi_c) ;
   
    // The spheres  
    Point center (ndim) ;
    center.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ; center.set(2) = 0 ; center.set(3) = 0 ;
    a_minus = aa*cosh(eta_minus)/sinh(eta_minus) ;
    domains[0] = new Domain_nucleus(0, ttype, r1, center, res) ;
    center.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    a_plus = aa*cosh(eta_plus)/sinh(eta_plus) ;
    domains[1] = new Domain_nucleus(1, ttype, r2, center, res) ;

    // Bispherical domains antitrigo order: 
    domains[2] = new Domain_bispheric_chi_first(2, ttype, aa, eta_minus, rr(0), chi_lim, res_bi) ;
    domains[3] = new Domain_bispheric_rect(3, ttype, aa, rr(0), eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[4] = new Domain_bispheric_eta_first(4, ttype, aa, rr(0), -eta_lim, eta_lim, res_bi) ;
    domains[5] = new Domain_bispheric_rect(5, ttype, aa, rr(0), eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[6] = new Domain_bispheric_chi_first(6, ttype, aa, eta_plus, rr(0), chi_lim, res_bi) ;
    
    // Shells     
    center.set(1) = 0 ;
    for (int i=0 ; i<nshells ; i++)
	domains[7+i] = new Domain_shell (7+i, ttype, rr(i), rr(i+1), center, res) ;

    // ZEC
    domains[7+nshells] = new Domain_compact (7+nshells, ttype,rr(nshells), center, res) ;
}

Space_bispheric::Space_bispheric (int ttype, double distance, double r1, double r2, int nn, const Array<double>& rr, const Array<int>& type_r, int nr) {

    // Verif :
    ndim = 3 ;
    
    ndom_minus = 1 ;
    ndom_plus = 1 ;
    nshells = nn ;
    
    nbr_domains = 8 + nshells ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
  
    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;
    
    // Bispheric :
    // Computation of aa
    Param par_a ;
    par_a.add_double(r1,0) ;
    par_a.add_double(r2,1) ;
    par_a.add_double(distance,2) ;
    double a_min = 0 ;
    double a_max = distance/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_a, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/r2) ;
    double eta_minus = -asinh(aa/r1) ;

    double chi_c = 2*atan(aa/rr(0)) ;
    double eta_c = log((1+rr(0)/aa)/(rr(0)/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rr(0), aa, chi_c) ;
   
    // The spheres  
    Point center (ndim) ;
    center.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ; center.set(2) = 0 ; center.set(3) = 0 ;
    a_minus = aa*cosh(eta_minus)/sinh(eta_minus) ;
    domains[0] = new Domain_nucleus(0, ttype, r1, center, res) ;
    center.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    a_plus = aa*cosh(eta_plus)/sinh(eta_plus) ;
    domains[1] = new Domain_nucleus(1, ttype, r2, center, res) ;

    // Bispherical domains antitrigo order: 
    domains[2] = new Domain_bispheric_chi_first(2, ttype, aa, eta_minus, rr(0), chi_lim, res_bi) ;
    domains[3] = new Domain_bispheric_rect(3, ttype, aa, rr(0), eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[4] = new Domain_bispheric_eta_first(4, ttype, aa, rr(0), -eta_lim, eta_lim, res_bi) ;
    domains[5] = new Domain_bispheric_rect(5, ttype, aa, rr(0), eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[6] = new Domain_bispheric_chi_first(6, ttype, aa, eta_plus, rr(0), chi_lim, res_bi) ;
    
    // Shells     
    center.set(1) = 0 ;
    for (int i=0 ; i<nshells ; i++) {
	switch (type_r(i)) {
	  case STD_TYPE: 
	    domains[7+i] = new Domain_shell (7+i,ttype, rr(i), rr(i+1), center, res) ;
	    break ;
	  case LOG_TYPE:
	    domains[7+i] = new Domain_shell_log (7+i,ttype, rr(i), rr(i+1), center, res) ;
	    break ;
	  case SURR_TYPE:
	    domains[7+i] = new Domain_shell_surr (7+i, ttype, rr(i), rr(i+1), center, res) ;
	    break ;
	  default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}
    }
    
    // ZEC
    domains[7+nshells] = new Domain_compact (7+nshells, ttype,rr(nshells), center, res) ;
}

Space_bispheric::Space_bispheric (int ttype, double distance, int nminus, const Array<double>& rminus, int nplus, const Array<double>& rplus, int nn, const Array<double>& rr, const Array<int>& type_r, int nr, bool withnuc) {

    // Verif :
    ndim = 3 ;
    
    ndom_minus = nminus ;
    ndom_plus = nplus ;
    nshells = nn ;
    
    nbr_domains = ndom_minus + ndom_plus + nshells + 6 ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
  
    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;
    
    // Bispheric :
    // Computation of aa
    double r1 = withnuc ? rminus(nminus-1) : rminus(nminus) ;
    double r2 = withnuc ? rplus(nplus-1) : rplus(nplus) ;
    
    if (fabs(r1-r2)>1e-12) {
      cerr << "Constructor of Space_bispheric not correct for different radii" << endl ;
      abort() ;
    }
    
    Param par_a ;
    par_a.add_double(r1,0) ;
    par_a.add_double(r2,1) ;
    par_a.add_double(distance,2) ;
    double a_min = 0 ;
    double a_max = distance/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_a, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/r2) ;
    double eta_minus = -asinh(aa/r1) ;

    double chi_c = 2*atan(aa/rr(0)) ;
    double eta_c = log((1+rr(0)/aa)/(rr(0)/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rr(0), aa, chi_c) ;
   
   // The spheres  
    int current = 0 ;    
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    a_minus = aa*cosh(eta_minus)/sinh(eta_minus) ;
     if (withnuc) {
    	domains[current] = new Domain_nucleus(current, ttype, rminus(0), center_minus, res) ;
    	current ++ ;

    for (int i=1 ; i<ndom_minus ; i++) {
	domains[current] = new Domain_shell (current, ttype, rminus(i-1), rminus(i), center_minus, res) ;
	 current ++ ;
    	}
    }
    else {
    for (int i=0 ; i<ndom_minus ; i++) {
	domains[current] = new Domain_shell (current, ttype, rminus(i), rminus(i+1), center_minus, res) ;
	current ++ ;
    }
    }

    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    a_plus = aa*cosh(eta_plus)/sinh(eta_plus) ;
    if (withnuc) {
  	domains[current] = new Domain_nucleus(current, ttype, rplus(0), center_plus, res) ;
  	current ++ ;
    
     for (int i=1 ; i<ndom_plus ; i++) {
	domains[current] = new Domain_shell (current, ttype, rplus(i-1), rplus(i), center_plus, res) ;
	current ++ ;
    }
    }
    else {
      for (int i=0 ; i<ndom_plus ; i++) {
	domains[current] = new Domain_shell (current, ttype, rplus(i), rplus(i+1), center_plus, res) ;
	current ++ ;
    }
    }
    

    // Bispherical domains antitrigo order: 
    domains[current] = new Domain_bispheric_chi_first(current, ttype, aa, eta_minus, rr(0), chi_lim, res_bi) ;
    domains[current+1] = new Domain_bispheric_rect(current+1, ttype, aa, rr(0), eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[current+2] = new Domain_bispheric_eta_first(current+2, ttype, aa, rr(0), -eta_lim, eta_lim, res_bi) ;
    domains[current+3] = new Domain_bispheric_rect(current+3,ttype, aa, rr(0), eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[current+4] = new Domain_bispheric_chi_first(current+4,ttype, aa, eta_plus, rr(0), chi_lim, res_bi) ;
    current += 5 ;
    // Shells     
    Point center(3) ;
    for (int i=0 ; i<nshells ; i++) {
	switch (type_r(i)) {
	  case STD_TYPE: 
	    domains[current] = new Domain_shell (current, ttype, rr(i), rr(i+1), center, res) ;
	    current ++ ;
	    break ;
	  case LOG_TYPE:
	    domains[current] = new Domain_shell_log (current, ttype, rr(i), rr(i+1), center, res) ;
	    current ++ ;
	    break ;
	  case SURR_TYPE:
	    domains[current] = new Domain_shell_surr (current, ttype, rr(i), rr(i+1), center, res) ;
            current ++ ;
	    break ;
	  default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}
    }
    
    // ZEC
    domains[current] = new Domain_compact (current, ttype,rr(nshells), center, res) ;
}

Space_bispheric::Space_bispheric (int ttype, double distance, int nminus, const Array<double>& rminus, const Array<int>& type_r_minus, int nplus, const Array<double>& rplus, const Array<int>& type_r_plus, int nn, const Array<double>& rr, const Array<int>& type_r, int nr, bool withnuc) {

    // Verif :
    ndim = 3 ;
    
    ndom_minus = nminus ;
    ndom_plus = nplus ;
    nshells = nn ;
    
    nbr_domains = ndom_minus + ndom_plus + nshells + 6  ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
  
    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;
    
    // Bispheric :
    // Computation of aa
    double r1 = rminus(nminus -1 ) ;
    double r2 = rplus(nplus-1) ;
   
    if (fabs(r1-r2)>1e-12) {
      cerr << "Constructor of Space_bispheric not correct for different radii" << endl ;
      abort() ;
    }
    
    Param par_a ;
    par_a.add_double(r1,0) ;
    par_a.add_double(r2,1) ;
    par_a.add_double(distance,2) ;
    double a_min = 0 ;
    double a_max = distance/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_a, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/r2) ;
    double eta_minus = -asinh(aa/r1) ;

    double chi_c = 2*atan(aa/rr(0)) ;
    double eta_c = log((1+rr(0)/aa)/(rr(0)/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rr(0), aa, chi_c) ;
   
    // The spheres  
    int current = 0 ;    
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    a_minus = aa*cosh(eta_minus)/sinh(eta_minus) ;
    if (withnuc) {
    	domains[current] = new Domain_nucleus(current, ttype, rminus(0), center_minus, res) ;
    	current ++ ;

    for (int i=1 ; i<ndom_minus ; i++) {
	switch (type_r_minus(i)) {
	  case STD_TYPE: 
	    domains[current] = new Domain_shell (current, ttype, rminus(i-1), rminus(i), center_minus, res) ;
	    current ++ ;
	    break ;
	  case LOG_TYPE:
	    domains[current] = new Domain_shell_log (current, ttype, rminus(i-1), rminus(i), center_minus, res) ;
	    current ++ ;
	    break ;
	  case SURR_TYPE:
	    domains[current] = new Domain_shell_surr (current, ttype, rminus(i-1), rminus(i), center_minus, res) ;
            current ++ ;
	    break ;
	  default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}
    }
    }
    else {
    for (int i=0 ; i<ndom_minus ; i++) {
	switch (type_r_minus(i)) {
	  case STD_TYPE: 
	    domains[current] = new Domain_shell (current, ttype, rminus(i), rminus(i+1), center_minus, res) ;
	    current ++ ;
	    break ;
	  case LOG_TYPE:
	    domains[current] = new Domain_shell_log (current, ttype, rminus(i), rminus(i+1), center_minus, res) ;
	    current ++ ;
	    break ;
	  case SURR_TYPE:
	    domains[current] = new Domain_shell_surr (current, ttype, rminus(i), rminus(i+1), center_minus, res) ;
            current ++ ;
	    break ;
	  default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}
    }
    }

    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    a_plus = aa*cosh(eta_plus)/sinh(eta_plus) ;
    if (withnuc) {
  	domains[current] = new Domain_nucleus(current, ttype, rplus(0), center_plus, res) ;
  	current ++ ;
    
     for (int i=1 ; i<ndom_plus ; i++) {
	switch (type_r_minus(i)) {
	  case STD_TYPE: 
	    domains[current] = new Domain_shell (current, ttype, rplus(i-1), rplus(i), center_plus, res) ;
	    current ++ ;
	    break ;
	  case LOG_TYPE:
	    domains[current] = new Domain_shell_log (current, ttype, rplus(i-1), rplus(i), center_plus, res) ;
	    current ++ ;
	    break ;
	  case SURR_TYPE:
	    domains[current] = new Domain_shell_surr (current, ttype, rplus(i-1), rplus(i), center_plus, res) ;
            current ++ ;
	    break ;
	  default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}
    }
    }
    else {
      for (int i=0 ; i<ndom_plus ; i++) {
	switch (type_r_minus(i)) {
	  case STD_TYPE: 
	    domains[current] = new Domain_shell (current, ttype, rplus(i), rplus(i+1), center_plus, res) ;
	    current ++ ;
	    break ;
	  case LOG_TYPE:
	    domains[current] = new Domain_shell_log (current, ttype, rplus(i), rplus(i+1), center_plus, res) ;
	    current ++ ;
	    break ;
	  case SURR_TYPE:
	    domains[current] = new Domain_shell_surr (current, ttype, rplus(i), rplus(i+1), center_plus, res) ;
            current ++ ;
	    break ;
	  default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}
    }
    }
    

    // Bispherical domains antitrigo order: 
    domains[current] = new Domain_bispheric_chi_first(current, ttype, aa, eta_minus, rr(0), chi_lim, res_bi) ;
    domains[current+1] = new Domain_bispheric_rect(current+1, ttype, aa, rr(0), eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[current+2] = new Domain_bispheric_eta_first(current+2, ttype, aa, rr(0), -eta_lim, eta_lim, res_bi) ;
    domains[current+3] = new Domain_bispheric_rect(current+3,ttype, aa, rr(0), eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[current+4] = new Domain_bispheric_chi_first(current+4,ttype, aa, eta_plus, rr(0), chi_lim, res_bi) ;
    current += 5 ;
    // Shells     
    Point center(3) ;
    for (int i=0 ; i<nshells ; i++) {
	switch (type_r(i)) {
	  case STD_TYPE: 
	    domains[current] = new Domain_shell (current, ttype, rr(i), rr(i+1), center, res) ;
	    current ++ ;
	    break ;
	  case LOG_TYPE:
	    domains[current] = new Domain_shell_log (current, ttype, rr(i), rr(i+1), center, res) ;
	    current ++ ;
	    break ;
	  case SURR_TYPE:
	    domains[current] = new Domain_shell_surr (current, ttype, rr(i), rr(i+1), center, res) ;
            current ++ ;
	    break ;
	  default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}
    }
    
    // ZEC
    domains[current] = new Domain_compact (current, ttype,rr(nshells), center, res) ;
}

Space_bispheric::Space_bispheric (int ttype, double distance, double rhor1, double rshell1, double rhor2, double rshell2, double rext, int nr) {

    // Verif :
    ndim = 3 ;
    
    ndom_minus = 1 ;
    ndom_plus = 1 ;
    nshells = 0 ;
    
    nbr_domains =8 ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
  
    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;
    
    // Bispheric :
    // Computation of aa
    double r1 = rshell1 ;
    double r2 = rshell2 ;
    
    Param par_a ;
    par_a.add_double(r1,0) ;
    par_a.add_double(r2,1) ;
    par_a.add_double(distance,2) ;
    double a_min = 0 ;
    double a_max = distance/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_a, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/r2) ;
    double eta_minus = -asinh(aa/r1) ;

    double chi_c = 2*atan(aa/rext) ;
    double eta_c = log((1+rext/aa)/(rext/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;
   
    // The spheres  
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    a_minus = aa*cosh(eta_minus)/sinh(eta_minus) ;

    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    a_plus = aa*cosh(eta_plus)/sinh(eta_plus) ;

    domains[0] = new Domain_shell(0, ttype, rhor1, rshell1, center_minus, res) ;
    domains[1] = new Domain_shell(1, ttype, rhor2, rshell2, center_plus, res) ;
   
    // Bispherical domains antitrigo order: 
    domains[2] = new Domain_bispheric_chi_first(2, ttype, aa, eta_minus, rext, chi_lim, res_bi) ;
    domains[3] = new Domain_bispheric_rect(3, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[4] = new Domain_bispheric_eta_first(4, ttype, aa, rext, -eta_lim, eta_lim, res_bi) ;
    domains[5] = new Domain_bispheric_rect(5,ttype, aa, rext, eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[6] = new Domain_bispheric_chi_first(6,ttype, aa, eta_plus, rext, chi_lim, res_bi) ;
    
    // Shells     
    Point center(3) ;
    domains[7] = new Domain_compact (7, ttype, rext, center, res) ;
}


Space_bispheric::Space_bispheric (int ttype, double distance, double rhor1, double rshell1, double rhor2, double rshell2, double rext, Dim_array** resol) {

    // Verif :
    ndim = 3 ;
    
    ndom_minus = 1 ;
    ndom_plus = 1 ;
    nshells = 0 ;
    
    nbr_domains =8 ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
  
    // Bispheric :
    // Computation of aa
    double r1 = rshell1 ;
    double r2 = rshell2 ;
    
    Param par_a ;
    par_a.add_double(r1,0) ;
    par_a.add_double(r2,1) ;
    par_a.add_double(distance,2) ;
    double a_min = 0 ;
    double a_max = distance/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_a, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/r2) ;
    double eta_minus = -asinh(aa/r1) ;

    double chi_c = 2*atan(aa/rext) ;
    double eta_c = log((1+rext/aa)/(rext/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;
   
    // The spheres  
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    a_minus = aa*cosh(eta_minus)/sinh(eta_minus) ;

    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    a_plus = aa*cosh(eta_plus)/sinh(eta_plus) ;

    domains[0] = new Domain_shell(0, ttype, rhor1, rshell1, center_minus, *resol[0]) ;
    domains[1] = new Domain_shell(1, ttype, rhor2, rshell2, center_plus, *resol[1]) ;
   
    // Bispherical domains antitrigo order: 
    domains[2] = new Domain_bispheric_chi_first(2, ttype, aa, eta_minus, rext, chi_lim, *resol[2]) ;
    domains[3] = new Domain_bispheric_rect(3, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, *resol[3]) ;
    domains[4] = new Domain_bispheric_eta_first(4, ttype, aa, rext, -eta_lim, eta_lim, *resol[4]) ;
    domains[5] = new Domain_bispheric_rect(5,ttype, aa, rext, eta_plus, eta_lim, chi_lim, *resol[5]) ;
    domains[6] = new Domain_bispheric_chi_first(6,ttype, aa, eta_plus, rext, chi_lim, *resol[6]) ;
    
    // Shells     
    Point center(3) ;
    domains[7] = new Domain_compact (7, ttype, rext, center, *resol[7]) ;
}

Space_bispheric::Space_bispheric (int ttype, double distance, double rhor1, double rshell1, double rhor2, double rshell2, int nn, const Array<double>& rr, int nr) {

    // Verif :
    ndim = 3 ;
    
    ndom_minus = 1 ;
    ndom_plus = 1 ;
    nshells = nn ;
    
    nbr_domains =8+nn ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
  
    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;
    
    // Bispheric :
    // Computation of aa
    double r1 = rshell1 ;
    double r2 = rshell2 ;

    double rext = rr(0) ;
    
    Param par_a ;
    par_a.add_double(r1,0) ;
    par_a.add_double(r2,1) ;
    par_a.add_double(distance,2) ;
    double a_min = 0 ;
    double a_max = distance/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_a, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/r2) ;
    double eta_minus = -asinh(aa/r1) ;

    double chi_c = 2*atan(aa/rext) ;
    double eta_c = log((1+rext/aa)/(rext/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;
   
    // The spheres  
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    a_minus = aa*cosh(eta_minus)/sinh(eta_minus) ;

    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    a_plus = aa*cosh(eta_plus)/sinh(eta_plus) ;

    domains[0] = new Domain_shell(0, ttype, rhor1, rshell1, center_minus, res) ;
    domains[1] = new Domain_shell(1, ttype, rhor2, rshell2, center_plus, res) ;
   
    // Bispherical domains antitrigo order: 
    domains[2] = new Domain_bispheric_chi_first(2, ttype, aa, eta_minus, rext, chi_lim, res_bi) ;
    domains[3] = new Domain_bispheric_rect(3, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[4] = new Domain_bispheric_eta_first(4, ttype, aa, rext, -eta_lim, eta_lim, res_bi) ;
    domains[5] = new Domain_bispheric_rect(5,ttype, aa, rext, eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[6] = new Domain_bispheric_chi_first(6,ttype, aa, eta_plus, rext, chi_lim, res_bi) ;
    
    // Shells       
    Point center(3) ;
    for (int i=0 ; i<nn ; i++)
	domains[7+i] = new Domain_shell (7+i, ttype, rr(i), rr(i+1), center, res) ;
    domains[7+nn] = new Domain_compact (7+nn, ttype, rr(nn), center, res) ;
}

Space_bispheric::Space_bispheric (FILE*fd, int type_shells, bool old) {
  
	fread_be (&nbr_domains, sizeof(int), 1, fd) ;
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;	
	fread_be (&a_minus, sizeof(double), 1, fd) ;
	fread_be (&a_plus, sizeof(double), 1, fd) ;
	
	if (!old) {
	  fread_be (&ndom_minus, sizeof(int), 1, fd) ;
	  fread_be (&ndom_plus, sizeof(int), 1, fd) ;
	  fread_be (&nshells, sizeof(int), 1, fd) ;
	}
	else {
	    ndom_minus = 1 ;
	    ndom_plus = 1 ;
	    nshells = nbr_domains-8 ;
	}

	// Check whether one has nucleii or not
	double nnuc = nbr_domains - 1 - nshells - ndom_minus - ndom_plus - 5;
	bool nucleus = (nnuc>=2) ? true : false ;
		
	domains = new Domain* [nbr_domains] ;

	int current = 0 ;
	if (nucleus) {
		domains[current] = new Domain_nucleus(current, fd) ;
		current ++ ;
	}

        for (int i=0 ; i<ndom_minus ; i++) {
	  domains[current] = new Domain_shell(current, fd) ;
	  current ++ ;
	}
	
	if (nucleus) {
		domains[current] = new Domain_nucleus(current, fd) ;
		current ++ ;
	}

	for (int i=0 ; i<ndom_plus ; i++) {
	  domains[current] = new Domain_shell(current, fd) ;
	  current ++ ;
 	}

    	// Bispherical domains antitrigo order: 
    	domains[current] = new Domain_bispheric_chi_first(current, fd) ;
	current ++ ;
    	domains[current] = new Domain_bispheric_rect(current, fd) ;
	current ++ ;
    	domains[current] = new Domain_bispheric_eta_first(current, fd) ;
	current ++ ;
    	domains[current] = new Domain_bispheric_rect(current, fd) ;
	current ++ ;
    	domains[current] = new Domain_bispheric_chi_first(current, fd) ;
	current ++ ;

	// Shells :
	for (int i=0 ; i<nshells ; i++)
	  switch (type_shells) {
	    case STD_TYPE :
	      domains[current] = new Domain_shell(current, fd) ;
		current ++ ;
	      break ;
	    case LOG_TYPE :
	      domains[current] = new Domain_shell_log(current, fd) ;
		current ++ ;
	      break ;
	    case SURR_TYPE :
	      domains[current] = new Domain_shell_surr(current, fd) ;
		current ++ ;
	      break ;
	   default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}

	// Compact
	domains[current] = new Domain_compact(current, fd) ;
	current ++ ;
}

Space_bispheric::Space_bispheric (FILE*fd, int type_minus, int type_plus, int type_shells, bool nucleus) {
  

	fread_be (&nbr_domains, sizeof(int), 1, fd) ;
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;	
	fread_be (&a_minus, sizeof(double), 1, fd) ;
	fread_be (&a_plus, sizeof(double), 1, fd) ;
	fread_be (&ndom_minus, sizeof(int), 1, fd) ;
	fread_be (&ndom_plus, sizeof(int), 1, fd) ;
	fread_be (&nshells, sizeof(int), 1, fd) ;
	
	domains = new Domain* [nbr_domains] ;

	int current = 0 ;
	if (nucleus) {
		domains[current] = new Domain_nucleus(current, fd) ;
		current ++ ;
	}

	int limminus = (nucleus) ? ndom_minus-1 : ndom_minus ;
	for (int i=0 ; i<limminus ; i++)
	  switch (type_minus) {
	    case STD_TYPE :
	      domains[current] = new Domain_shell(current, fd) ;
		current ++ ;
	      break ;
	    case LOG_TYPE :
	      domains[current] = new Domain_shell_log(current, fd) ;
		current ++ ;
	      break ;
	    case SURR_TYPE :
	      domains[current] = new Domain_shell_surr(current, fd) ;
		current ++ ;
	      break ;
	   default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}
	
	if (nucleus) {
		domains[current] = new Domain_nucleus(current, fd) ;
		current ++ ;
	}
	
	int limplus = (nucleus) ? ndom_plus-1 : ndom_plus ;
	for (int i=0 ; i<limplus ; i++)
	  switch (type_plus) {
	    case STD_TYPE :
	      domains[current] = new Domain_shell(current, fd) ;
		current ++ ;
	      break ;
	    case LOG_TYPE :
	      domains[current] = new Domain_shell_log(current, fd) ;
		current ++ ;
	      break ;
	    case SURR_TYPE :
	      domains[current] = new Domain_shell_surr(current, fd) ;
		current ++ ;
	      break ;
	   default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}

    	// Bispherical domains antitrigo order: 
    	domains[current] = new Domain_bispheric_chi_first(current, fd) ;
	current ++ ;
    	domains[current] = new Domain_bispheric_rect(current, fd) ;
	current ++ ;
    	domains[current] = new Domain_bispheric_eta_first(current, fd) ;
	current ++ ;
    	domains[current] = new Domain_bispheric_rect(current, fd) ;
	current ++ ;
    	domains[current] = new Domain_bispheric_chi_first(current, fd) ;
	current ++ ;

	// Shells :
	for (int i=0 ; i<nshells ; i++)
	  switch (type_shells) {
	    case STD_TYPE :
	      domains[current] = new Domain_shell(current, fd) ;
		current ++ ;
	      break ;
	    case LOG_TYPE :
	      domains[current] = new Domain_shell_log(current, fd) ;
		current ++ ;
	      break ;
	    case SURR_TYPE :
	      domains[current] = new Domain_shell_surr(current, fd) ;
		current ++ ;
	      break ;
	   default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}

	// Compact
	domains[current] = new Domain_compact(current, fd) ;
	current ++ ;
}

Space_bispheric::~Space_bispheric() 
{
}


void Space_bispheric::save (FILE* fd) const  {
	fwrite_be (&nbr_domains, sizeof(int), 1, fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	fwrite_be (&a_minus, sizeof(double), 1, fd) ;
	fwrite_be (&a_plus, sizeof(double), 1, fd) ;
	fwrite_be (&ndom_minus, sizeof(int), 1, fd) ;
	fwrite_be (&ndom_plus, sizeof(int), 1, fd) ;
	fwrite_be (&nshells, sizeof(int), 1, fd) ;
	for (int i=0 ; i<nbr_domains ; i++)
		domains[i]->save(fd) ;
}

Array<int> Space_bispheric::get_indices_matching_non_std(int dom, int bound) const {

	if (dom==ndom_minus-1) {
	  // First sphere ;
	  Array<int> res (2,2) ;
	  switch (bound) {
		case OUTER_BC :
		  res.set(0,0) = ndom_minus+ndom_plus ; // Matching with chi first
		  res.set(1,0) = INNER_BC ;
		  res.set(0,1) = ndom_minus+ndom_plus+1 ; // Matching with rect
		  res.set(1,1) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bispheric::get_indices_matching_non_std" << endl ;
		  abort() ;
	  }
	return res ;
	}
	
	if (dom== ndom_minus+ndom_plus-1) {
		// second sphere ;
		Array<int> res(2, 2) ;
		switch (bound) {
		case OUTER_BC :
		  res.set(0,0) = ndom_minus+ndom_plus+3 ; // Matching with rect
		  res.set(1,0) = INNER_BC ;
		  res.set(0,1) = ndom_minus+ndom_plus+4 ; // Matching with chi_first
		  res.set(1,1) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bispheric::get_indices_matching_non_std" << endl ;			abort() ;
		}
	return res ;
	}

	if (dom== ndom_minus+ndom_plus)	{
	  // first chi first :
	  Array<int> res(2,1) ;
	  switch (bound) {
	    case INNER_BC : 
	      res.set(0,0) = ndom_minus-1 ; // First sphere
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = ndom_minus+ndom_plus+5 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bispheric::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}
	
	if (dom==ndom_minus+ndom_plus+1)	{
	  // first rect :
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case INNER_BC : 
	      res.set(0,0) = ndom_minus-1 ; // First sphere
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = ndom_minus+ndom_plus+5 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bispheric::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}

	    
	if (dom==ndom_plus+ndom_minus+2) {
	  // eta first 
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case OUTER_BC :
	      res.set(0, 0) = ndom_minus+ndom_plus+5 ; // Compactified domain or first shell
	      res.set(1, 0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bispheric::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}	
	
	if (dom==ndom_minus+ndom_plus+3) {
	  // second rect :
	  Array<int> res(2, 1) ;
	  switch (bound) {
		case INNER_BC : 
		  res.set(0, 0) = ndom_minus+ndom_plus-1 ; // Second sphere
		  res.set(1, 0) = OUTER_BC ;
		  break ;
		case OUTER_BC :
		  res.set(0, 0) = ndom_minus+ndom_plus+5 ; // Compactified domain or first shell
		  res.set(1, 0) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bispheric::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}

	if (dom==ndom_minus+ndom_plus+4) {
	  // second chi first :
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case INNER_BC : 
	      res.set(0,0) = ndom_minus+ndom_plus-1 ; // second nucleus
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = ndom_minus+ndom_plus+5 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
		cerr << "Bad bound in Space_bispheric::get_indices_matching_non_std" << endl ;
		abort() ;
	  }
	return res ;
	}
	
	if (dom==ndom_minus+ndom_plus+5) {
	    // compactified domain or first shell :
	    Array<int> res(2,5) ;
	    switch (bound) {
		case INNER_BC :
		  res.set(0,0) = ndom_minus+ndom_plus ; // first chi first
		  res.set(0,1) = ndom_minus+ndom_plus+1 ; // first rect
		  res.set(0,2) = ndom_minus+ndom_plus+2 ; // eta first
		  res.set(0,3) = ndom_minus+ndom_plus+3 ; // second rect
		  res.set(0,4) = ndom_minus+ndom_plus+4 ; // second chi first
		  // Outer BC for all :
		  for (int i=0 ; i<5 ; i++)
		      res.set(1,i) = OUTER_BC ;
		  break ;
	  default :
	      cerr << "Bad bound in Space_bispheric::get_indices_matching_non_std" << endl ;
	      abort() ;
	}
	return res ;
	}
			
	cerr << "Bad domain in Space_bispheric::get_indices_matching_non_std" << endl ;
	abort() ;
}


double Space_bispheric::int_inf (const Scalar& so) const {

	const Domain_compact* p_cmp = dynamic_cast <const Domain_compact*> (domains[nbr_domains-1]) ;
	if (p_cmp==0x0) {
		cerr << "No compactified domain in Space_bispheric::int_inf" << endl ;
		abort() ;
	}
	return p_cmp->integ (so(nbr_domains-1), OUTER_BC) ;
}


double Space_bispheric::int_sphere_one (const Scalar& so) const {

	if (ndom_minus==1) {
	  return (domains[ndom_minus+ndom_plus]->integ(so(ndom_minus+ndom_plus), INNER_BC) + domains[ndom_minus+ndom_plus+1]->integ(so(ndom_minus+ndom_plus+1), INNER_BC)) ;
	}
	else {
	  return domains[1]->integ(so(1), INNER_BC) ;
	}
}

double Space_bispheric::int_sphere_two (const Scalar& so) const {

	if (ndom_plus==1) {
	  return (domains[ndom_minus+ndom_plus+3]->integ(so(ndom_minus+ndom_plus+3), INNER_BC) + domains[ndom_minus+ndom_plus+4]->integ(so(ndom_minus+ndom_plus+4), INNER_BC)) ;
	}
	else {
	    return domains[ndom_minus+1]->integ(so(1), INNER_BC) ;
	}
}

}
