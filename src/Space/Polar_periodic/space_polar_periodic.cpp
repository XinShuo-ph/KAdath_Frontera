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
#include "polar_periodic.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "system_of_eqs.hpp"

namespace Kadath {
Space_polar_periodic::Space_polar_periodic(int ttype, double omega, const Dim_array& res, const Array<double>& bounds) {

    // Verif :
    assert (bounds.get_ndim()== 1) ;

    ndim = 3 ;
    
    nbr_domains = bounds.get_size(0) ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
    // Nucleus
    domains[0] = new Domain_polar_periodic_nucleus(0, ttype, bounds(0), omega, res) ;
    for (int i=1 ; i<nbr_domains ; i++)
       domains[i] = new Domain_polar_periodic_shell(i, ttype, bounds(i-1), bounds(i), omega, res) ;
  
}

Space_polar_periodic::Space_polar_periodic(FILE* fd) {
	fread_be (&nbr_domains, sizeof(int), 1, fd) ;
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;
	domains = new Domain* [nbr_domains] ;
	//nucleus :
	domains[0] = new Domain_polar_periodic_nucleus(0, fd) ;
	//Shells :
	for (int i=1 ; i<nbr_domains ; i++)
		domains[i] = new Domain_polar_periodic_shell(i, fd) ;
}

Space_polar_periodic::~Space_polar_periodic() {
}

void Space_polar_periodic::save (FILE* fd) const  {
	fwrite_be (&nbr_domains, sizeof(int), 1, fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;	
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	for (int i=0 ; i<nbr_domains ; i++)
		domains[i]->save(fd) ;
}

/*
int Space_polar_periodic::nbr_unknowns_from_variable_domains () const {
  
    return domains[0]->nbr_unknowns_from_adapted() ;
}

void Space_polar_periodic::affecte_coef_to_variable_domains (int& pos, int cc, Array<int>& zedoms) const {
  int old_pos = pos ;
  bool indic = false ;
  bool firstindic ;
  for (int d=0 ; d<nbr_domains ; d++) {
	pos = old_pos ;
	domains[d]->affecte_coef(pos, cc, indic) ;
	if (d==0)
		firstindic = indic ;
	assert(indic==firstindic) ;
	if (indic) 
		zedoms.set(d) = d ;
	else
		zedoms.set(d) = -1 ;  
	}
}


void Space_polar_periodic::xx_to_ders_variable_domains (const Array<double>& xx, int& conte) const {
  
  int old_conte = conte ;
  for (int d=0 ; d<nbr_domains ; d++) {
	  conte = old_conte ;
	  domains[d]->xx_to_ders_from_adapted(xx, conte) ;
  }
}


void Space_polar_periodic::xx_to_vars_variable_domains (System_of_eqs* sys, const Array<double>& xx, int& pos) const  {
  
  // First get the correction :
  int old_pos = pos ;
  double cor_ome ;
  domains[0]->xx_to_vars_from_adapted(cor_ome, xx, pos) ;

  // Now update the variables :
  for (int i=0 ; i<sys->nvar ; i++) { 
    for (int n=0 ; n<sys->var[i]->get_n_comp() ; n++) {
	Scalar res(*this) ;
        for (int d=0 ; d<nbr_domains ; d++)
		domains[d]->update_variable(cor_ome, *sys->var[i]->cmp[n], res) ;

	for (int d=0 ; d<nbr_domains ; d++)
		sys->var[i]->cmp[n]->set_domain(d) = res(d) ;
    }
  }
  
  // Now the constants :
  for (int i=0 ; i<sys->ncst ; i++) 
    if (sys->cst[i*(sys->dom_max-sys->dom_min+1)]->get_type_data()==TERM_T)
      for (int n=0 ; n<sys->cst[i*(sys->dom_max-sys->dom_min+1)]->val_t->get_n_comp() ; n++) {
	
	Scalar so (*this) ;
	so = 0 ;
	for (int d=0 ; d<nbr_domains ; d++)
	  so.set_domain(d) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1)+d]->val_t->cmp[n])(d) ;
	
	Scalar res(*this) ;
	res = 0 ;
	for (int d=0 ; d<nbr_domains ; d++)
		domains[d]->update_constante(cor_ome, so, res) ;
	for (int d=0 ; d<nbr_domains ; d++) 
		sys->cst[i*(sys->dom_max-sys->dom_min+1)+1+d]->val_t->cmp[n]->set_domain(d) = res(d) ;
      }
      
      // Update the mapping :
      for (int d=0 ; d<nbr_domains ; d++)
	      domains[d]->update_mapping(cor_ome) ;
}
*/

double Space_polar_periodic::get_omega() const {
	  const Domain_polar_periodic_nucleus* polarperiodicnuc = dynamic_cast<const Domain_polar_periodic_nucleus*>(domains[0]) ;
	return polarperiodicnuc->get_ome() ;
     
}
}
