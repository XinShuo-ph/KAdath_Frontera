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

Space_bispheric::Space_bispheric (int ttype, double distance, int nminus, const Array<double>& rminus, int nplus, const Array<double>& rplus, int nn, const Array<double>& rr, const Array<int>& type_r, int nr) {

    // Verif :
    ndim = 3 ;
    
    ndom_minus = nminus ;
    ndom_plus = nplus ;
    nshells = nn ;
    
    nbr_domains = ndom_minus + ndom_plus + nshells + 6;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
  
    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;
    
    // Bispheric :
    // Computation of aa
    double r1 = rminus(nminus-1) ;
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
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    a_minus = aa*cosh(eta_minus)/sinh(eta_minus) ;
    domains[0] = new Domain_nucleus(0, ttype, rminus(0), center_minus, res) ;
    
    for (int i=1 ; i<ndom_minus ; i++)
      domains[i] = new Domain_shell(i, ttype, rminus(i-1), rminus(i), center_minus, res) ;
      
    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    a_plus = aa*cosh(eta_plus)/sinh(eta_plus) ;
    domains[ndom_minus] = new Domain_nucleus(ndom_minus, ttype, rplus(0), center_plus, res) ;
  
 
    for (int i=1 ; i<ndom_plus ; i++)
       domains[i+ndom_minus] = new Domain_shell(i+ndom_minus, ttype, rplus(i-1), rplus(i), center_plus, res) ;
    
    // Bispherical domains antitrigo order: 
    domains[ndom_minus+ndom_plus] = new Domain_bispheric_chi_first(ndom_minus+ndom_plus, ttype, aa, eta_minus, rr(0), chi_lim, res_bi) ;
    domains[ndom_minus+ndom_plus+1] = new Domain_bispheric_rect(ndom_minus+ndom_plus+1, ttype, aa, rr(0), eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[ndom_minus+ndom_plus+2] = new Domain_bispheric_eta_first(ndom_minus+ndom_plus+2, ttype, aa, rr(0), -eta_lim, eta_lim, res_bi) ;
    domains[ndom_minus+ndom_plus+3] = new Domain_bispheric_rect(ndom_minus+ndom_plus+3,ttype, aa, rr(0), eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[ndom_minus+ndom_plus+4] = new Domain_bispheric_chi_first(ndom_minus+ndom_plus+4,ttype, aa, eta_plus, rr(0), chi_lim, res_bi) ;
    
    // Shells     
    Point center(3) ;
    for (int i=0 ; i<nshells ; i++) {
	switch (type_r(i)) {
	  case STD_TYPE: 
	    domains[ndom_plus+ndom_minus+5+i] = new Domain_shell (ndom_plus+ndom_minus+5+i, ttype, rr(i), rr(i+1), center, res) ;
	    break ;
	  case LOG_TYPE:
	    domains[ndom_plus+ndom_minus+5+i] = new Domain_shell_log (ndom_plus+ndom_minus+5+i, ttype, rr(i), rr(i+1), center, res) ;
	    break ;
	  case SURR_TYPE:
	    domains[ndom_plus+ndom_minus+5+i] = new Domain_shell_surr (ndom_plus+ndom_minus+5+i, ttype, rr(i), rr(i+1), center, res) ;
	    break ;
	  default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}
    }
    
    // ZEC
    domains[nbr_domains-1] = new Domain_compact (nbr_domains-1, ttype,rr(nshells), center, res) ;
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
	
	domains = new Domain* [nbr_domains] ;
	
	domains[0] = new Domain_nucleus(0, fd) ;

        for (int i=1 ; i<ndom_minus ; i++)
	  domains[i] = new Domain_shell(i, fd) ;
	
	domains[ndom_minus] = new Domain_nucleus(ndom_minus, fd) ;

	for (int i=ndom_minus+1 ; i<ndom_minus+ndom_plus ; i++)
	  domains[i] = new Domain_shell(i, fd) ;
	

    	// Bispherical domains antitrigo order: 
    	domains[ndom_minus+ndom_plus] = new Domain_bispheric_chi_first(ndom_minus+ndom_plus, fd) ;
    	domains[ndom_minus+ndom_plus+1] = new Domain_bispheric_rect(ndom_minus+ndom_plus+1, fd) ;
    	domains[ndom_minus+ndom_plus+2] = new Domain_bispheric_eta_first(ndom_minus+ndom_plus+2, fd) ;
    	domains[ndom_minus+ndom_plus+3] = new Domain_bispheric_rect(ndom_minus+ndom_plus+3, fd) ;
    	domains[ndom_minus+ndom_plus+4] = new Domain_bispheric_chi_first(ndom_minus+ndom_plus+4, fd) ;

	// Shells :
	for (int i=0 ; i<nshells ; i++)
	  switch (type_shells) {
	    case STD_TYPE :
	      domains[i+ndom_minus+ndom_plus+5] = new Domain_shell(i+ndom_minus+ndom_plus+5, fd) ;
	      break ;
	    case LOG_TYPE :
	      domains[i+ndom_minus+ndom_plus+5] = new Domain_shell_log(i+ndom_minus+ndom_plus+5, fd) ;
	      break ;
	    case SURR_TYPE :
	      domains[i+ndom_minus+ndom_plus+5] = new Domain_shell_surr(i+ndom_minus+ndom_plus+5, fd) ;
	      break ;
	   default:
	     cerr << "Unknown type of shells in Space_bishperic constructor" << endl ;
	     abort() ;
	}

	// Compact
	domains[nbr_domains-1] = new Domain_compact(nbr_domains-1, fd) ;
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

	if (dom== ndom_plus+ndom_plus)	{
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
	
	if (dom==ndom_plus+ndom_plus+1)	{
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
