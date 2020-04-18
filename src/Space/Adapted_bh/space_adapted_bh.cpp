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
#include "adapted_bh.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"
namespace Kadath {
Space_adapted_bh::Space_adapted_bh (int ttype, const Point& center, const Dim_array& res, const Array<double>& bounds) {

    // Verif :
    assert (bounds.get_ndim()== 1) ;

    ndim = 3 ;
    
    nbr_domains = bounds.get_size(0) ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
     
    domains[0] = new Domain_shell_inner_adapted (*this, 1, ttype, bounds(0), bounds(1), center, res) ;
    for (int i=1 ; i<nbr_domains-1 ; i++)
       domains[i] = new Domain_shell(i, ttype, bounds(i), bounds(i+1), center, res) ;
  
   domains[nbr_domains-1] = new Domain_compact (nbr_domains-1, ttype, bounds(nbr_domains-1), center, res) ;
   
   const Domain_shell_inner_adapted* pshell_inner = dynamic_cast<const Domain_shell_inner_adapted*> (domains[0]) ;
   pshell_inner->vars_to_terms() ;
}


Space_adapted_bh::Space_adapted_bh (FILE* fd) {
	fread_be (&nbr_domains, sizeof(int), 1, fd) ;
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;
	domains = new Domain* [nbr_domains] ;
		
	domains[0] = new Domain_shell_inner_adapted (*this, 0, fd) ;
	//Shells :
	for (int i=1 ; i<nbr_domains-1 ; i++)
		domains[i] = new Domain_shell(i, fd) ;
	// Compactified
	domains[nbr_domains-1] = new Domain_compact(nbr_domains-1, fd) ;  
	
	const Domain_shell_inner_adapted* pshell_inner = dynamic_cast<const Domain_shell_inner_adapted*> (domains[0]) ;
	pshell_inner->vars_to_terms() ;
}

Space_adapted_bh::~Space_adapted_bh() {
  const Domain_shell_inner_adapted* pshell_inner = dynamic_cast<const Domain_shell_inner_adapted*> (domains[0]) ;
  pshell_inner->del_deriv() ;
}

void Space_adapted_bh::save (FILE* fd) const  {
	fwrite_be (&nbr_domains, sizeof(int), 1, fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;	
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	for (int i=0 ; i<nbr_domains ; i++)
		domains[i]->save(fd) ;
}

int Space_adapted_bh::nbr_unknowns_from_variable_domains () const {
    return domains[0]->nbr_unknowns_from_adapted() ;
}

void Space_adapted_bh::affecte_coef_to_variable_domains (int& pos, int cc, Array<int>& zedoms) const {
  
  bool found = false ;
  domains[0]->affecte_coef(pos, cc, found) ;
  
  if (found)
    zedoms.set(0) = 0 ;
  else
    zedoms.set(0) = -1 ;
}

void Space_adapted_bh::xx_to_ders_variable_domains (const Array<double>& xx, int& conte) const {
  
  domains[0]->xx_to_ders_from_adapted(xx, conte) ;
}

void Space_adapted_bh::xx_to_vars_variable_domains (System_of_eqs* sys, const Array<double>& xx, int& pos) const  {
  
  // First get the corrections :
  Val_domain cor_inner (domains[0]) ;
  domains[0]->xx_to_vars_from_adapted(cor_inner, xx, pos) ;
     
  // Now update the variables :
  for (int i=0 ; i<sys->nvar ; i++) { 
    for (int n=0 ; n<sys->var[i]->get_n_comp() ; n++) {
	Scalar res(*this) ;
	domains[0]->update_variable(cor_inner, *sys->var[i]->cmp[n], res) ;
	
	sys->var[i]->cmp[n]->set_domain(0) = res(0) ;
    }
  }
  
  // Now the constants :
  for (int i=0 ; i<sys->ncst ; i++) 
    if (sys->cst[i*(sys->dom_max-sys->dom_min+1)]->get_type_data()==TERM_T)
      for (int n=0 ; n<sys->cst[i*(sys->dom_max-sys->dom_min+1)]->val_t->get_n_comp() ; n++) {
	
	Scalar so (*this) ;
	so = 0 ;
	so.set_domain(0) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1)]->val_t->cmp[n])(0) ;
	
	Scalar res(*this) ;
	res = 0 ;
	domains[0]->update_constante(cor_inner, so, res) ;
	
	sys->cst[i*(sys->dom_max-sys->dom_min+1)]->val_t->cmp[n]->set_domain(0) = res(0) ;
      }
      
      // Update the mapping :
      domains[0]->update_mapping(cor_inner) ;
}

}
