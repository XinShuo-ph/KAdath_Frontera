/*
    Copyright 2022 Philippe Grandclement

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
#include "adapted_homothetic.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"
namespace Kadath {
Space_spheric_adapted_homothetic::Space_spheric_adapted_homothetic (int ttype, const Point& center, const Dim_array& res, const Array<double>& bounds) {

    // Verif :
    assert (bounds.get_ndim()== 1) ;

    ndim = 3 ;
    
    nbr_domains = bounds.get_size(0)+1 ;
    assert (nbr_domains>=4) ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
    
    // Nucleus
    domains[0] = new Domain_nucleus (0, ttype, bounds(0), center, res) ; 
    domains[1] = new Domain_shell_outer_homothetic (*this, 1, ttype, bounds(0), bounds(1), center, res) ;    
    domains[2] = new Domain_shell_inner_homothetic (*this, 2, ttype, bounds(1), bounds(2), center, res) ;
    for (int i=3 ; i<nbr_domains-1 ; i++)
       domains[i] = new Domain_shell(i, ttype, bounds(i-1), bounds(i), center, res) ;
  
   domains[nbr_domains-1] = new Domain_compact (nbr_domains-1, ttype, bounds(nbr_domains-2), center, res) ;
   
   const Domain_shell_outer_homothetic* pshell_outer = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[1]) ;
   pshell_outer->vars_to_terms() ;
   const Domain_shell_inner_homothetic* pshell_inner = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[2]) ;
   pshell_inner->vars_to_terms() ;
}


Space_spheric_adapted_homothetic::Space_spheric_adapted_homothetic (FILE* fd) {
	fread_be (&nbr_domains, sizeof(int), 1, fd) ;
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;
	domains = new Domain* [nbr_domains] ;
	//nucleus :
	domains[0] = new Domain_nucleus (0, fd) ;
	domains[1] = new Domain_shell_outer_homothetic (*this, 1, fd) ;	
	domains[2] = new Domain_shell_inner_homothetic (*this, 2, fd) ;
	//Shells :
	for (int i=3 ; i<nbr_domains-1 ; i++)
		domains[i] = new Domain_shell(i, fd) ;
	// Compactified
	domains[nbr_domains-1] = new Domain_compact(nbr_domains-1, fd) ;  
	
	const Domain_shell_outer_homothetic* pshell_outer = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[1]) ;
	pshell_outer->vars_to_terms() ;
	const Domain_shell_inner_homothetic* pshell_inner = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[2]) ;
	pshell_inner->vars_to_terms() ;
}

Space_spheric_adapted_homothetic::~Space_spheric_adapted_homothetic() {
    Domain_shell_outer_homothetic* pshell_outer = dynamic_cast<Domain_shell_outer_homothetic*> (domains[1]) ;
    pshell_outer->del_deriv() ;
    Domain_shell_inner_homothetic* pshell_inner = dynamic_cast<Domain_shell_inner_homothetic*> (domains[2]) ;
    pshell_inner->del_deriv() ;
}

void Space_spheric_adapted_homothetic::save (FILE* fd) const  {
	fwrite_be (&nbr_domains, sizeof(int), 1, fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;	
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	for (int i=0 ; i<nbr_domains ; i++)
		domains[i]->save(fd) ;
}

int Space_spheric_adapted_homothetic::nbr_unknowns_from_variable_domains () const {
  
    return domains[1]->nbr_unknowns_from_adapted() ;
}

void Space_spheric_adapted_homothetic::affecte_coef_to_variable_domains (int& pos, int cc, Array<int>& zedoms) const {
  
  bool found_outer = false ;
  int old_pos = pos ;
  domains[1]->affecte_coef(pos, cc, found_outer) ;
  pos = old_pos ;
  bool found_inner = false ;
  domains[2]->affecte_coef(pos, cc, found_inner) ;
  assert (found_outer == found_inner) ;
  if (found_outer) {
    zedoms.set(0) = 1 ;
    zedoms.set(1) = 2 ;
  }
  else {
    zedoms.set(0) = -1 ;
    zedoms.set(1) = -1 ;
  }
}

void Space_spheric_adapted_homothetic::xx_to_ders_variable_domains (const Array<double>& xx, int& conte) const {
  
  int old_conte = conte ;
  domains[1]->xx_to_ders_from_adapted(xx, conte) ;
  conte = old_conte ;
  domains[2]->xx_to_ders_from_adapted(xx, conte) ;
}

void Space_spheric_adapted_homothetic::xx_to_vars_variable_domains (System_of_eqs* sys, const Array<double>& xx, int& pos) const  {
  
  // First get the corrections :
  int old_pos = pos ;
  Val_domain cor_outer (domains[1]) ;
  domains[1]->xx_to_vars_from_adapted(cor_outer, xx, pos) ;
  pos = old_pos ;
  Val_domain cor_inner (domains[2]) ;
  domains[2]->xx_to_vars_from_adapted(cor_inner, xx, pos) ;
      
  // Now update the variables :
  for (int i=0 ; i<sys->nvar ; i++) { 
    for (int n=0 ; n<sys->var[i]->get_n_comp() ; n++) {
	Scalar res(*this) ;
	domains[1]->update_variable(cor_outer, *sys->var[i]->cmp[n], res) ;
	domains[2]->update_variable(cor_inner, *sys->var[i]->cmp[n], res) ;
	
	sys->var[i]->cmp[n]->set_domain(1) = res(1) ;
	sys->var[i]->cmp[n]->set_domain(2) = res(2) ; 
	
    }
  }
  
  // Now the constants :
  for (int i=0 ; i<sys->ncst ; i++) 
    if (sys->cst[i*(sys->dom_max-sys->dom_min+1)]->get_type_data()==TERM_T)
      for (int n=0 ; n<sys->cst[i*(sys->dom_max-sys->dom_min+1)]->val_t->get_n_comp() ; n++) {
	
	Scalar so (*this) ;
	so = 0 ;
	for (int d=1 ; d<=2 ; d++)
	  so.set_domain(d) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1)+d-sys->dom_min]->val_t->cmp[n])(d) ;
	  
	Scalar res(*this) ;
	res = 0 ;
	domains[1]->update_constante(cor_outer, so, res) ;
	domains[2]->update_constante(cor_inner, so, res) ;
	
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+1-sys->dom_min]->val_t->cmp[n]->set_domain(1) = res(1) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+2-sys->dom_min]->val_t->cmp[n]->set_domain(2) = res(2) ;
	
      }
      
      // Update the mapping :
      domains[1]->update_mapping(cor_outer) ;
      domains[2]->update_mapping(cor_inner) ;
}


Array<int> Space_spheric_adapted_homothetic::get_indices_matching_non_std(int dom, int bound) const {

	assert ((dom>=0) && (dom<nbr_domains)) ;
	Array<int> res (2, 1) ;
	switch (bound) {
		case OUTER_BC : 
			res.set(0,0) = dom+1 ;
			res.set(1,0) = INNER_BC ; 
			break ;
		case INNER_BC :
			res.set(0,0) = dom-1 ;
			res.set(1,0) = OUTER_BC ;
			break ;
		default :
			cerr << "Unknown boundary in " << endl ;
			cerr << *this << endl ;
			abort() ;
		}
	return res ;
}
}
