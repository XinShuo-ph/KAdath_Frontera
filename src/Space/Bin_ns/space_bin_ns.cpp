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
#include "bin_ns.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"
namespace Kadath {
double eta_lim_chi (double chi, double rext, double a, double eta_c) ;
double chi_lim_eta (double chi, double rext, double a, double chi_c) ;

double zerosec(double (*f)(double, const Param&), const Param& parf, 
    double x1, double x2, double precis, int nitermax, int& niter) ;



double func_abns (double aa, const Param& par) {
	double r1 = par.get_double(0) ;
	double r2 = par.get_double(1) ;
	double d = par.get_double(2) ;
	return (sqrt(aa*aa+r1*r1)+sqrt(aa*aa+r2*r2)-d) ;
}

Space_bin_ns::Space_bin_ns (int ttype, double dist, double rinstar1, double rstar1, double routstar1, 
			double rinstar2, double rstar2, double routstar2, double rext, int nr) {

    ndim = 3 ;
    
    nshells = 0 ;
    nbr_domains = 12 ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
    
    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;
  
     // Bispheric :
    // Computation of aa
    Param par_a ;
    par_a.add_double(routstar1,0) ;
    par_a.add_double(routstar2,1) ;
    par_a.add_double(dist,2) ;
    double a_min = 0 ;
    double a_max = dist/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_abns, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/routstar2) ;
    double eta_minus = -asinh(aa/routstar1) ;
    
    double chi_c = 2*atan(aa/rext) ;
    double eta_c = log((1+rext/aa)/(rext/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;
    
    // First NS
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    domains[0] = new Domain_nucleus (0, ttype, rinstar1, center_minus, res) ; 
    domains[1] = new Domain_shell_outer_adapted (*this, 1, ttype, rinstar1, rstar1, center_minus, res) ;    
    domains[2] = new Domain_shell_inner_adapted (*this, 2, ttype, rstar1, routstar1, center_minus, res) ;
   
    // Second NS 
    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    domains[3] = new Domain_nucleus (3, ttype, rinstar2, center_plus, res) ; 
    domains[4] = new Domain_shell_outer_adapted (*this, 4, ttype, rinstar2, rstar2, center_plus, res) ;    
    domains[5] = new Domain_shell_inner_adapted (*this, 5, ttype, rstar2, routstar2, center_plus, res) ;
    
    // Bispheric part
    domains[6] = new Domain_bispheric_chi_first(6, ttype, aa, eta_minus, rext, chi_lim, res_bi) ;
    domains[7] = new Domain_bispheric_rect(7, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[8] = new Domain_bispheric_eta_first(8, ttype, aa, rext, -eta_lim, eta_lim, res_bi) ;
    domains[9] = new Domain_bispheric_rect(9, ttype, aa, rext, eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[10] = new Domain_bispheric_chi_first(10, ttype, aa, eta_plus, rext, chi_lim, res_bi) ;
    
    // Compactified domain
    Point center(3) ; 
    domains[11] = new Domain_compact (11, ttype, rext, center, res) ;
    
    const Domain_shell_outer_adapted* pouter_1 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[1]) ;
    pouter_1->vars_to_terms() ;
    pouter_1->update() ;
    const Domain_shell_inner_adapted* pinner_1 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[2]) ;
    pinner_1->vars_to_terms() ;
    pinner_1->update() ;
    const Domain_shell_outer_adapted* pouter_2 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[4]) ;
    pouter_2->vars_to_terms() ;
    pouter_2->update() ;
    const Domain_shell_inner_adapted* pinner_2 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[5]) ;
    pinner_2->vars_to_terms() ;
    pinner_2->update() ;
}

Space_bin_ns::Space_bin_ns (int ttype, double dist, double rinstar1, double rstar1, double routstar1, 
			double rinstar2, double rstar2, double routstar2, double rext, double rshell, int nr) {

    ndim = 3 ;
    
    nshells = 1 ;
    nbr_domains = 13 ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
    
    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;
  
     // Bispheric :
    // Computation of aa
    Param par_a ;
    par_a.add_double(routstar1,0) ;
    par_a.add_double(routstar2,1) ;
    par_a.add_double(dist,2) ;
    double a_min = 0 ;
    double a_max = dist/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_abns, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/routstar2) ;
    double eta_minus = -asinh(aa/routstar1) ;
    
    double chi_c = 2*atan(aa/rext) ;
    double eta_c = log((1+rext/aa)/(rext/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;
    
    // First NS
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    domains[0] = new Domain_nucleus (0, ttype, rinstar1, center_minus, res) ; 
    domains[1] = new Domain_shell_outer_adapted (*this, 1, ttype, rinstar1, rstar1, center_minus, res) ;    
    domains[2] = new Domain_shell_inner_adapted (*this, 2, ttype, rstar1, routstar1, center_minus, res) ;
   
    // Second NS 
    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    domains[3] = new Domain_nucleus (3, ttype, rinstar2, center_plus, res) ; 
    domains[4] = new Domain_shell_outer_adapted (*this, 4, ttype, rinstar2, rstar2, center_plus, res) ;    
    domains[5] = new Domain_shell_inner_adapted (*this, 5, ttype, rstar2, routstar2, center_plus, res) ;
    
    // Bispheric part
    domains[6] = new Domain_bispheric_chi_first(6, ttype, aa, eta_minus, rext, chi_lim, res_bi) ;
    domains[7] = new Domain_bispheric_rect(7, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[8] = new Domain_bispheric_eta_first(8, ttype, aa, rext, -eta_lim, eta_lim, res_bi) ;
    domains[9] = new Domain_bispheric_rect(9, ttype, aa, rext, eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[10] = new Domain_bispheric_chi_first(10, ttype, aa, eta_plus, rext, chi_lim, res_bi) ;
    
    // Shell in 1/R 
    Point center(3) ; 
    domains[11] = new Domain_shell (11, ttype, rext, rshell,  center, res) ;
  
    // Compactified domain
    domains[12] = new Domain_compact (12, ttype, rshell,  center, res) ;
 
    const Domain_shell_outer_adapted* pouter_1 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[1]) ;
    pouter_1->vars_to_terms() ;
    pouter_1->update() ;
    const Domain_shell_inner_adapted* pinner_1 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[2]) ;
    pinner_1->vars_to_terms() ;
    pinner_1->update() ;
    const Domain_shell_outer_adapted* pouter_2 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[4]) ;
    pouter_2->vars_to_terms() ;
    pouter_2->update() ;
    const Domain_shell_inner_adapted* pinner_2 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[5]) ;
    pinner_2->vars_to_terms() ;
    pinner_2->update() ;
}

Space_bin_ns::Space_bin_ns (FILE* fd) {
	fread_be (&nbr_domains, sizeof(int), 1, fd) ;
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;

	nshells = nbr_domains-12 ;

	domains = new Domain* [nbr_domains] ;
	//First NS :
	domains[0] = new Domain_nucleus (0, fd) ;
	domains[1] = new Domain_shell_outer_adapted (*this, 1, fd) ;	
	domains[2] = new Domain_shell_inner_adapted (*this, 2, fd) ;
	
	//second NS :
	domains[3] = new Domain_nucleus (3, fd) ;
	domains[4] = new Domain_shell_outer_adapted (*this, 4, fd) ;	
	domains[5] = new Domain_shell_inner_adapted (*this, 5, fd) ;
	
	// Bispheric
	domains[6] = new Domain_bispheric_chi_first(6, fd) ;
    	domains[7] = new Domain_bispheric_rect(7, fd) ;
    	domains[8] = new Domain_bispheric_eta_first(8, fd) ;
    	domains[9] = new Domain_bispheric_rect(9, fd) ;
    	domains[10] = new Domain_bispheric_chi_first(10, fd) ;
	
	for (int i=0 ; i<nshells ; i++)
		domains[11+i] = new Domain_shell (11+i, fd) ;

	// Compactified
	domains[11+nshells] = new Domain_compact(11+nshells, fd) ;  
	
    const Domain_shell_outer_adapted* pouter_1 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[1]) ;
    pouter_1->vars_to_terms() ;
    pouter_1->update() ;
    const Domain_shell_inner_adapted* pinner_1 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[2]) ;
    pinner_1->vars_to_terms() ;
    pinner_1->update() ;
    const Domain_shell_outer_adapted* pouter_2 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[4]) ;
    pouter_2->vars_to_terms() ;
    pouter_2->update() ;
    const Domain_shell_inner_adapted* pinner_2 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[5]) ;
    pinner_2->vars_to_terms() ;
    pinner_2->update() ;
    
}

Space_bin_ns::~Space_bin_ns() {
    Domain_shell_outer_adapted* pouter_1 = dynamic_cast<Domain_shell_outer_adapted*> (domains[1]) ;
    pouter_1->del_deriv() ;
    Domain_shell_inner_adapted* pinner_1 = dynamic_cast<Domain_shell_inner_adapted*> (domains[2]) ;
    pinner_1->del_deriv() ;
    Domain_shell_outer_adapted* pouter_2 = dynamic_cast<Domain_shell_outer_adapted*> (domains[4]) ;
    pouter_2->del_deriv() ;
    Domain_shell_inner_adapted* pinner_2 = dynamic_cast<Domain_shell_inner_adapted*> (domains[5]) ;
    pinner_2->del_deriv() ;
}

void Space_bin_ns::save (FILE* fd) const  {
	fwrite_be (&nbr_domains, sizeof(int), 1, fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;	
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	for (int i=0 ; i<nbr_domains ; i++)
		domains[i]->save(fd) ;
}

int Space_bin_ns::nbr_unknowns_from_variable_domains () const {
  
    return (domains[1]->nbr_unknowns_from_adapted() + domains[4]->nbr_unknowns_from_adapted()) ;
}

void Space_bin_ns::affecte_coef_to_variable_domains (int& pos, int cc, Array<int>& zedoms) const {
  
  // In star 1 ?
  bool found_outer_1 = false ;
  int old_pos = pos ;
  domains[1]->affecte_coef(pos, cc, found_outer_1) ;
  pos = old_pos ;
  bool found_inner_1 = false ;
  domains[2]->affecte_coef(pos, cc, found_inner_1) ;
  assert (found_outer_1 == found_inner_1) ;
  if (found_outer_1) {
    zedoms.set(0) = 1 ;
    zedoms.set(1) = 2 ;
  }
  else {
    zedoms.set(0) = -1 ;
    zedoms.set(1) = -1 ;
  }
  
  // In star 2 ?
    bool found_outer_2 = false ;
    old_pos = pos ;
    domains[4]->affecte_coef(pos, cc, found_outer_2) ;
    pos = old_pos ;
    bool found_inner_2 = false ;
    domains[5]->affecte_coef(pos, cc, found_inner_2) ;
    assert (found_outer_2 == found_inner_2) ;
    if (found_outer_2) {
      zedoms.set(0) = 4 ;
      zedoms.set(1) = 5 ;
    }
}

void Space_bin_ns::xx_to_ders_variable_domains (const Array<double>& xx, int& conte) const {
  
  // Star 1 
  int old_conte = conte ;
  domains[1]->xx_to_ders_from_adapted(xx, conte) ;
  conte = old_conte ;
  domains[2]->xx_to_ders_from_adapted(xx, conte) ;
  // Star 2
  old_conte = conte ;
  domains[4]->xx_to_ders_from_adapted(xx, conte) ;
  conte = old_conte ;
  domains[5]->xx_to_ders_from_adapted(xx, conte) ;
}

void Space_bin_ns::xx_to_vars_variable_domains (System_of_eqs* sys, const Array<double>& xx, int& pos) const  {
  

  // First get the corrections :
  // Star 1 
  int old_pos = pos ;
  Val_domain cor_outer_1 (domains[1]) ;
  domains[1]->xx_to_vars_from_adapted(cor_outer_1, xx, pos) ;
  pos = old_pos ;
  Val_domain cor_inner_1 (domains[2]) ;
  domains[2]->xx_to_vars_from_adapted(cor_inner_1, xx, pos) ;
      
  // Star 2
  old_pos = pos ;
  Val_domain cor_outer_2 (domains[4]) ;
  domains[4]->xx_to_vars_from_adapted(cor_outer_2, xx, pos) ;
  pos = old_pos ;
  Val_domain cor_inner_2 (domains[5]) ;
  domains[5]->xx_to_vars_from_adapted(cor_inner_2, xx, pos) ;
    
  // Now update the variables :
  for (int i=0 ; i<sys->nvar ; i++) { 
    for (int n=0 ; n<sys->var[i]->get_n_comp() ; n++) {
	Scalar res(*this) ;
	domains[1]->update_variable(cor_outer_1, *sys->var[i]->cmp[n], res) ;
	domains[2]->update_variable(cor_inner_1, *sys->var[i]->cmp[n], res) ;
	domains[4]->update_variable(cor_outer_2, *sys->var[i]->cmp[n], res) ;
	domains[5]->update_variable(cor_inner_2, *sys->var[i]->cmp[n], res) ;
	
	
	sys->var[i]->cmp[n]->set_domain(1) = res(1) ;
	sys->var[i]->cmp[n]->set_domain(2) = res(2) ; 
	sys->var[i]->cmp[n]->set_domain(4) = res(4) ;
	sys->var[i]->cmp[n]->set_domain(5) = res(5) ; 
    }
  }
  
  // Now the constants :
  for (int i=0 ; i<sys->ncst ; i++) 
    if (sys->cst[i*(sys->dom_max-sys->dom_min+1)]->get_type_data()==TERM_T)
      for (int n=0 ; n<sys->cst[i*(sys->dom_max-sys->dom_min+1)]->val_t->get_n_comp() ; n++) {
	
	Scalar so (*this) ;
	so = 0 ;
	for (int d=1 ; d<=5 ; d++)
	  if (d!=3)
	  so.set_domain(d) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1)+d]->val_t->cmp[n])(d) ;
	
	Scalar res(*this) ;
	res = 0 ;
	domains[1]->update_constante(cor_outer_1, so, res) ;
	domains[2]->update_constante(cor_inner_1, so, res) ;
	domains[4]->update_constante(cor_outer_2, so, res) ;
	domains[5]->update_constante(cor_inner_2, so, res) ;
	
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+1]->val_t->cmp[n]->set_domain(1) = res(1) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+2]->val_t->cmp[n]->set_domain(2) = res(2) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+4]->val_t->cmp[n]->set_domain(4) = res(4) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+5]->val_t->cmp[n]->set_domain(5) = res(5) ;
      }
      
      // Update the mapping :
      domains[1]->update_mapping(cor_outer_1) ;
      domains[2]->update_mapping(cor_inner_1) ; 
      domains[4]->update_mapping(cor_outer_2) ;
      domains[5]->update_mapping(cor_inner_2) ;
}


Array<int> Space_bin_ns::get_indices_matching_non_std(int dom, int bound) const {

	if (dom==2) {
	  // First sphere ;
	  Array<int> res (2,2) ;
	  switch (bound) {
		case OUTER_BC :
		  res.set(0,0) = 6 ; // Matching with chi first
		  res.set(1,0) = INNER_BC ;
		  res.set(0,1) = 7 ; // Matching with rect
		  res.set(1,1) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;
		  abort() ;
	  }
	return res ;
	}
	
	if (dom== 5) {
		// second sphere ;
		Array<int> res(2, 2) ;
		switch (bound) {
		case OUTER_BC :
		  res.set(0,0) = 9 ; // Matching with rect
		  res.set(1,0) = INNER_BC ;
		  res.set(0,1) = 10 ; // Matching with chi_first
		  res.set(1,1) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;			abort() ;
		}
	return res ;
	}

	if (dom== 6)	{
	  // first chi first :
	  Array<int> res(2,1) ;
	  switch (bound) {
	    case INNER_BC : 
	      res.set(0,0) = 2 ; // First sphere
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = 11 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}
	
	if (dom==7)	{
	  // first rect :
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case INNER_BC : 
	      res.set(0,0) = 2 ; // First sphere
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = 11 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}

	    
	if (dom==8) {
	  // eta first 
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case OUTER_BC :
	      res.set(0, 0) = 11 ; // Compactified domain or first shell
	      res.set(1, 0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}	
	
	if (dom==9) {
	  // second rect :
	  Array<int> res(2, 1) ;
	  switch (bound) {
		case INNER_BC : 
		  res.set(0, 0) = 5 ; // Second sphere
		  res.set(1, 0) = OUTER_BC ;
		  break ;
		case OUTER_BC :
		  res.set(0, 0) = 11 ; // Compactified domain or first shell
		  res.set(1, 0) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}

	if (dom==10) {
	  // second chi first :
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case INNER_BC : 
	      res.set(0,0) = 5 ; // second nucleus
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = 11 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
		cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;
		abort() ;
	  }
	return res ;
	}
	
	if (dom==11) {
	    // compactified domain or first shell :
	    Array<int> res(2,5) ;
	    switch (bound) {
		case INNER_BC :
		  res.set(0,0) = 6 ; // first chi first
		  res.set(0,1) = 7 ; // first rect
		  res.set(0,2) = 8 ; // eta first
		  res.set(0,3) = 9 ; // second rect
		  res.set(0,4) = 10 ; // second chi first
		  // Outer BC for all :
		  for (int i=0 ; i<5 ; i++)
		      res.set(1,i) = OUTER_BC ;
		  break ;
	  default :
	      cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;
	      abort() ;
	}
	return res ;
	}
			
	cerr << "Bad domain in Space_bispheric::get_indices_matching_non_std" << endl ;
	abort() ;
}

}
