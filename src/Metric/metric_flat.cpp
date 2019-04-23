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

#include "space.hpp"
#include "system_of_eqs.hpp"
#include "tensor.hpp"
#include "metric.hpp"
#include "term_eq.hpp"
#include "scalar.hpp"
#include "name_tools.hpp"
#include "metric_tensor.hpp"
namespace Kadath {
Metric_flat::Metric_flat (const Space& sp, const Base_tensor& bb) : Metric(sp), basis(bb) {
}

Metric_flat::Metric_flat (const Metric_flat& so) : Metric (so), basis(so.basis) {
}

Metric_flat::~Metric_flat() {
}

void Metric_flat::update() {
	// Nothing to do everything is constant
}

void Metric_flat::update(int) {
	// Nothing to do everything is constant
}

void Metric_flat::compute_cov (int dd) const {
	Metric_tensor res (espace, COV, basis) ;

	for (int i=1 ; i<=espace.get_ndim()  ; i++)
	  for (int j=i ; j<=espace.get_ndim() ; j++)
	    res.set(i,j).set_domain(dd) = (i==j) ? 1. : 0 ;
	res.std_base() ;
	
	if (p_met_cov[dd]==0x0)
		p_met_cov[dd] = new Term_eq(dd, res) ;
	else
		*p_met_cov[dd] = Term_eq (dd, res) ;
	p_met_cov[dd]->set_der_zero() ;
}

int Metric_flat::give_type(int dd) const {
  return basis.get_basis(dd) ;
}

void Metric_flat::compute_con (int dd) const {

	Metric_tensor res (espace, CON, basis) ;
	for (int i=1 ; i<=espace.get_ndim() ; i++)
	  for (int j=i ; j<=espace.get_ndim() ; j++)
	    res.set(i,j).set_domain(dd) = (i==j) ? 1. : 0 ;
	res.std_base() ;
	if (p_met_con[dd]==0x0)
		p_met_con[dd] = new Term_eq(dd, res) ;
	else
		*p_met_con[dd] = Term_eq(dd, res) ;
	p_met_con[dd]->set_der_zero() ;
}

void Metric_flat::compute_christo (int) const {
  cerr << "Computation of Christo not explicit for Metric_flat" << endl ;
  abort() ;
}

void Metric_flat::compute_ricci_tensor (int dd) const {

	Metric_tensor res (espace, COV, basis) ;
	switch (basis.get_basis(dd)) {
		case CARTESIAN_BASIS :
			res = 0 ;
			break ;
		case SPHERICAL_BASIS :
			res = 0 ;
			break ;
		case MTZ_BASIS :
			for (int i=1 ; i<=3 ; i++)
				for (int j=i ; j<=3 ; j++)
					res.set(i,j).annule_hard() ;
			res.set(2,2).set_domain(dd) = -2/pow(espace.get_domain(dd)->get_radius(), 2) ;
			res.set(3,3).set_domain(dd) = res(2,2)(dd) ;
			break ;
		default:
			cerr << "Unknown tensorial basis in Metric_flat::compute_ricci_tensor" << endl ;
			abort() ;
	}

	res.std_base() ;
	if (p_ricci_tensor[dd]==0x0)
		p_ricci_tensor[dd] = new Term_eq(dd, res) ;
	else
		*p_ricci_tensor[dd] = Term_eq(dd, res) ;

	p_ricci_tensor[dd]->set_der_zero() ;
}

void Metric_flat::compute_ricci_scalar (int dd) const {

	Scalar res (espace) ;
	switch (basis.get_basis(dd)) {
		case CARTESIAN_BASIS :
			res = 0 ;
			break ;
		case SPHERICAL_BASIS :
			res = 0 ;
			break ;
		case MTZ_BASIS :
			res.set_domain(dd)= -4/pow(espace.get_domain(dd)->get_radius(), 2) ;
			break ;
		default:
			cerr << "Unknown tensorial basis in Metric_flat::compute_ricci_tensor" << endl ;
			abort() ;
	}

	res.std_base() ;
	if (p_ricci_scalar[dd]==0x0)
		p_ricci_scalar[dd] = new Term_eq(dd, res) ;
	else
		*p_ricci_scalar[dd] = Term_eq(dd, res) ;

	p_ricci_scalar[dd]->set_der_zero() ;
}

void Metric_flat::manipulate_ind (Term_eq& so, int ind) const {
	// Just change the type of the indice !
	so.set_val_t()->set_index_type (ind) *= -1 ;
	if (so.set_der_t() !=0x0)
		so.set_der_t()->set_index_type (ind) *= -1 ;
}


Term_eq Metric_flat::derive_partial_cart (int type_der, char ind_der, const Term_eq& so) const {
  
	int dom = so.get_dom() ;
	
	// Affecte names or not ?
	bool donames = ((so.val_t->is_name_affected()) || (so.val_t->get_valence()==0)) ? true : false ;
	
	// Computation of flat gradient :
	Term_eq auxi (espace.get_domain(dom)->partial_cart (so)) ;
	
	int val_res = auxi.val_t->get_valence() ;
	
	if (donames) {
	// Set the names of the indices :
	auxi.val_t->set_name_affected() ;
	auxi.val_t->set_name_ind(0, ind_der) ;
	for (int i=1 ; i<val_res ; i++)
		auxi.val_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;

	if (auxi.der_t!=0x0) {
	  auxi.der_t->set_name_affected() ;
	  auxi.der_t->set_name_ind(0, ind_der) ;
	  for (int i=1 ; i<val_res ; i++)
		auxi.der_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;
	}
	}
	
	// Manipulate if Contravariant version
	if (type_der==CON) 
		manipulate_ind (auxi, 0) ;
	
	bool need_sum = false ;
	if (donames)
	  for (int i=1 ; i<val_res ; i++)
		if (ind_der== so.val_t->get_name_ind()[i-1])
			need_sum = true ;
	
	if (!need_sum)
		return auxi ;
	else {
	  if (auxi.der_t==0x0)
	    return Term_eq (dom, auxi.val_t->do_summation_one_dom(dom)) ;
	  else
	    return Term_eq (dom, auxi.val_t->do_summation_one_dom(dom), auxi.der_t->do_summation_one_dom(dom)) ;
	}
}

Term_eq Metric_flat::derive_partial_spher (int type_der, char ind_der, const Term_eq& so) const {
  
	int dom = so.get_dom() ;
	bool donames = ((so.val_t->is_name_affected()) || (so.val_t->get_valence()==0)) ? true : false ;
	
	// Computation of flat gradient :
	Term_eq auxi (espace.get_domain(dom)->partial_spher (so)) ;
	
	int val_res = auxi.val_t->get_valence() ;
	
	if (donames) {
	// Set the names of the indices :
	auxi.val_t->set_name_affected() ;
	auxi.val_t->set_name_ind(0, ind_der) ;
	for (int i=1 ; i<val_res ; i++)
		auxi.val_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;

	if (auxi.der_t!=0x0) {
	  auxi.der_t->set_name_affected() ;
	  auxi.der_t->set_name_ind(0, ind_der) ;
	  for (int i=1 ; i<val_res ; i++)
		auxi.der_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;
	}
	}
	
	// Manipulate if Contravariant version
	if (type_der==CON) 
		manipulate_ind (auxi, 0) ;
	
	bool need_sum = false ;
	if (donames)
	for (int i=1 ; i<val_res ; i++)
		if (ind_der== so.val_t->get_name_ind()[i-1])
			need_sum = true ;
	
	if (!need_sum)
		return auxi ;
	else {
	  if (auxi.der_t==0x0)
	    return Term_eq (dom, auxi.val_t->do_summation_one_dom(dom)) ;
	  else
	    return Term_eq (dom, auxi.val_t->do_summation_one_dom(dom), auxi.der_t->do_summation_one_dom(dom)) ;
	}
}

Term_eq Metric_flat::derive_partial_mtz (int type_der, char ind_der, const Term_eq& so) const {
 
	int dom = so.get_dom() ;
	bool donames = ((so.val_t->is_name_affected()) || (so.val_t->get_valence()==0)) ? true : false ;
	
	// Computation of flat gradient :
	Term_eq auxi (espace.get_domain(dom)->partial_mtz (so)) ;
	
	int val_res = auxi.val_t->get_valence() ;
	
	if (donames) {
	// Set the names of the indices :
	auxi.val_t->set_name_affected() ;
	auxi.val_t->set_name_ind(0, ind_der) ;
	for (int i=1 ; i<val_res ; i++)
		auxi.val_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;

	if (auxi.der_t!=0x0) {
	  auxi.der_t->set_name_affected() ;
	  auxi.der_t->set_name_ind(0, ind_der) ;
	  for (int i=1 ; i<val_res ; i++)
		auxi.der_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;
	}
	}
	
	// Manipulate if Contravariant version
	if (type_der==CON) 
		manipulate_ind (auxi, 0) ;
	
	bool need_sum = false ;
	if (donames)
	for (int i=1 ; i<val_res ; i++)
		if (ind_der== so.val_t->get_name_ind()[i-1])
			need_sum = true ;
	
	if (!need_sum)
		return auxi ;
	else {
	  if (auxi.der_t==0x0)
	    return Term_eq (dom, auxi.val_t->do_summation_one_dom(dom)) ;
	  else
	    return Term_eq (dom, auxi.val_t->do_summation_one_dom(dom), auxi.der_t->do_summation_one_dom(dom)) ;
	}
}



Term_eq Metric_flat::derive_partial (int type_der, char ind_der, const Term_eq& so) const {
  
	int dom = so.get_dom() ;
      
	if (p_met_con[dom]==0x0)
		compute_con(dom) ;
	if (p_met_cov[dom]==0x0)
		compute_cov(dom) ;

	// so must be tensor :
	if (so.get_type_data()!=TERM_T) {
		cerr << "Metric_flat::derive partial only defined for tensor data" << endl ;
		abort() ;
	}

	switch (basis.get_basis(dom)) {
	  case CARTESIAN_BASIS :
	      return derive_partial_cart (type_der, ind_der, so) ;
	  case SPHERICAL_BASIS :
	      return derive_partial_spher (type_der, ind_der, so) ;
	  case MTZ_BASIS :
	      return derive_partial_mtz (type_der, ind_der, so) ;
	  default:
	      cerr << "Unknown tensorial basis in Metric_flat::derive_partial" << endl ;
	      abort() ;
	}
}

Term_eq Metric_flat::derive_cart (int type_der, char ind_der, const Term_eq& so) const {
  return derive_partial (type_der, ind_der, so) ;
}

Term_eq Metric_flat::derive_spher (int type_der, char ind_der, const Term_eq& so) const {

	// Computation of flat gradient :
	Term_eq part_der (derive_partial(type_der, ind_der, so)) ;
	
	int dom = so.get_dom() ;
	bool donames = ((so.val_t->is_name_affected()) || (so.val_t->get_valence()==0)) ? true : false ;
	
	Term_eq auxi (espace.get_domain(dom)->connection_spher (so)) ;   
	int val_res = auxi.val_t->get_valence() ;
	
	if (donames) {
	// Set the names of the indices :
	auxi.val_t->set_name_affected() ;
	auxi.val_t->set_name_ind(0, ind_der) ;
	for (int i=1 ; i<val_res ; i++)
		auxi.val_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;

	if (auxi.der_t!=0x0) {
	  auxi.der_t->set_name_affected() ;
	  auxi.der_t->set_name_ind(0, ind_der) ;
	  for (int i=1 ; i<val_res ; i++)
		auxi.der_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;
	}
	}
	
	// Manipulate if Contravariant version
	if (type_der==CON) 
		manipulate_ind (auxi, 0) ;
	
	
	bool need_sum = false ;
	if (donames) 
	for (int i=1 ; i<val_res ; i++)
		if (ind_der== so.val_t->get_name_ind()[i-1])
			need_sum = true ;
	
	if (!need_sum) {
	    return (part_der + auxi) ;
	}
	else {
	  if (auxi.der_t==0x0) {
	   
	   
	    return (part_der +  Term_eq (dom, auxi.val_t->do_summation_one_dom(dom))) ;
	  } 
	  else
	    return (part_der + Term_eq (dom, auxi.val_t->do_summation_one_dom(dom), auxi.der_t->do_summation_one_dom(dom))) ;
	}
}

Term_eq Metric_flat::derive_mtz (int type_der, char ind_der, const Term_eq& so) const {


	// Computation of flat gradient :
	Term_eq part_der (derive_partial(type_der, ind_der, so)) ;

	int dom = so.get_dom() ;
	bool donames = ((so.val_t->is_name_affected()) || (so.val_t->get_valence()==0)) ? true : false ;
	
	Term_eq auxi (espace.get_domain(dom)->connection_mtz (so)) ;   
	int val_res = auxi.val_t->get_valence() ;
	
	if (donames) {
	// Set the names of the indices :
	auxi.val_t->set_name_affected() ;
	auxi.val_t->set_name_ind(0, ind_der) ;
	for (int i=1 ; i<val_res ; i++)
		auxi.val_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;

	if (auxi.der_t!=0x0) {
	  auxi.der_t->set_name_affected() ;
	  auxi.der_t->set_name_ind(0, ind_der) ;
	  for (int i=1 ; i<val_res ; i++)
		auxi.der_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;
	}
	}
	
	// Manipulate if Contravariant version
	if (type_der==CON) 
		manipulate_ind (auxi, 0) ;
	
	
	bool need_sum = false ;
	if (donames) 
	for (int i=1 ; i<val_res ; i++)
		if (ind_der== so.val_t->get_name_ind()[i-1])
			need_sum = true ;
	
	if (!need_sum) {
	    return (part_der + auxi) ;
	}
	else {
	  if (auxi.der_t==0x0) {
	   
	   
	    return (part_der +  Term_eq (dom, auxi.val_t->do_summation_one_dom(dom))) ;
	  } 
	  else
	    return (part_der + Term_eq (dom, auxi.val_t->do_summation_one_dom(dom), auxi.der_t->do_summation_one_dom(dom))) ;
	}
}

Term_eq Metric_flat::derive (int type_der, char ind_der, const Term_eq& so) const {

	int dom = so.get_dom() ;

	if (p_met_con[dom]==0x0)
		compute_con(dom) ;
	if (p_met_cov[dom]==0x0)
		compute_cov(dom) ;

	// so must be tensor :
	if (so.get_type_data()!=TERM_T) {
		cerr << "Metric_flat::derive only defined for tensor data" << endl ;
		abort() ;
	}

      switch (basis.get_basis(dom)) {
	  case CARTESIAN_BASIS :
	      return derive_cart (type_der, ind_der, so) ;
	  case SPHERICAL_BASIS :
	      return derive_spher (type_der, ind_der, so) ;
	  case MTZ_BASIS :
	      return derive_mtz (type_der, ind_der, so) ;
	  default:
	      cerr << "Unknown tensorial basis in Metric_flat::derive" << endl ;
	      abort() ;
	}
	
}
Term_eq Metric_flat::derive_with_other_cart (int type_der, char ind_der, const Term_eq& so, const Metric* manipulator) const {
        int dom = so.get_dom() ;
	// Call the domain version
	return so.val_t->get_space().get_domain(dom)->derive_flat_cart (type_der, ind_der, so, manipulator) ;
}

Term_eq Metric_flat::derive_with_other_spher (int type_der, char ind_der, const Term_eq& so, const Metric* manipulator) const {
        int dom = so.get_dom() ;
	
	// Call the domain version
	return so.val_t->get_space().get_domain(dom)->derive_flat_spher (type_der, ind_der, so, manipulator) ;
}

Term_eq Metric_flat::derive_with_other_mtz (int type_der, char ind_der, const Term_eq& so, const Metric* manipulator) const {
        int dom = so.get_dom() ;
	
	// Call the domain version
	return so.val_t->get_space().get_domain(dom)->derive_flat_mtz (type_der, ind_der, so, manipulator) ;
}

Term_eq Metric_flat::derive_with_other (int type_der, char ind_der, const Term_eq& so, const Metric* manipulator) const {

	int dom = so.get_dom() ;

	if (p_met_con[dom]==0x0)
		compute_con(dom) ;
	if (p_met_cov[dom]==0x0)
		compute_cov(dom) ;

	// so must be tensor :
	if (so.get_type_data()!=TERM_T) {
		cerr << "Metric_flat::derive_with_other only defined for tensor data" << endl ;
		abort() ;
	}

      switch (basis.get_basis(dom)) {
	  case CARTESIAN_BASIS :
	      return derive_with_other_cart (type_der, ind_der, so, manipulator) ;
	  case SPHERICAL_BASIS :
	      return derive_with_other_spher (type_der, ind_der, so, manipulator) ; 
	case MTZ_BASIS :
	      return derive_with_other_mtz (type_der, ind_der, so, manipulator) ;
	  default:
	      cerr << "Unknown tensorial basis in Metric_flat::derive_with_other" << endl ;
	      abort() ;
	}
	
}

void Metric_flat::set_system (System_of_eqs& ss, const char* name) {
	
	syst = &ss ;
	if (syst->met!=0x0) {
		cerr << "Metric already set for the system" << endl ;
		abort() ;
	}
	
	ss.met = this ;
	ss.name_met = new char[LMAX] ;
	trim_spaces (ss.name_met, name) ;
}}
