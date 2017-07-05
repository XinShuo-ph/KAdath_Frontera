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
#include "tensor.hpp"
#include "metric.hpp"
#include "term_eq.hpp"
#include "scalar.hpp"
#include "system_of_eqs.hpp"
namespace Kadath {
Metric::Metric (const Space& sp) : espace (sp), syst(0x0), type_tensor(0) {
	p_met_cov = new Term_eq* [espace.get_nbr_domains()] ;
	p_met_con = new Term_eq* [espace.get_nbr_domains()] ;
	p_christo = new Term_eq* [espace.get_nbr_domains()] ;	
	p_riemann = new Term_eq* [espace.get_nbr_domains()] ;
	p_ricci_tensor = new Term_eq* [espace.get_nbr_domains()] ;
	p_ricci_scalar = new Term_eq* [espace.get_nbr_domains()] ;
	p_dirac = new Term_eq* [espace.get_nbr_domains()] ;
	p_det_cov = new Term_eq* [espace.get_nbr_domains()] ;
	for (int i=0 ; i<espace.get_nbr_domains() ; i++) {
		p_met_con[i] = 0x0 ;
		p_met_cov[i] = 0x0 ;
		p_christo[i] = 0x0 ;	
		p_riemann[i] = 0x0 ;
		p_ricci_tensor[i] = 0x0 ;
		p_ricci_scalar[i] = 0x0 ;
		p_dirac[i] = 0x0 ;
		p_det_cov[i] = 0x0 ;
	}
}

Metric::Metric (const Metric& so) : espace (so.espace), syst(so.syst), type_tensor(so.type_tensor) {

	p_met_cov = new Term_eq* [espace.get_nbr_domains()] ;
	p_met_con = new Term_eq* [espace.get_nbr_domains()] ;
	p_christo = new Term_eq* [espace.get_nbr_domains()] ;
	p_riemann = new Term_eq* [espace.get_nbr_domains()] ;
	p_ricci_tensor = new Term_eq* [espace.get_nbr_domains()] ;
	p_ricci_scalar = new Term_eq* [espace.get_nbr_domains()] ;
	p_dirac = new Term_eq* [espace.get_nbr_domains()] ;
	p_det_cov = new Term_eq* [espace.get_nbr_domains()] ;
	for (int i=0 ; i<espace.get_nbr_domains() ; i++) {
		p_met_cov[i] = (so.p_met_cov[i]==0x0) ? 0x0 : new Term_eq(*so.p_met_cov[i]) ;
		p_met_con[i] = (so.p_met_con[i]==0x0) ? 0x0 : new Term_eq(*so.p_met_con[i]) ;
		p_christo[i] = (so.p_christo[i]==0x0) ? 0x0 : new Term_eq(*so.p_christo[i]) ;
		p_riemann[i] = (so.p_riemann[i]==0x0) ? 0x0 : new Term_eq(*so.p_riemann[i]) ;
		p_ricci_tensor[i] = (so.p_ricci_tensor[i]==0x0) ? 0x0 : new Term_eq(*so.p_ricci_tensor[i]) ;
		p_ricci_scalar[i] = (so.p_ricci_scalar[i]==0x0) ? 0x0 : new Term_eq(*so.p_ricci_scalar[i]) ;
		p_dirac[i] = (so.p_dirac[i]==0x0) ? 0x0 : new Term_eq(*so.p_dirac[i]) ;
		p_det_cov[i] = (so.p_det_cov[i]==0x0) ? 0x0 : new Term_eq(*so.p_det_cov[i]) ;
	}
}

Metric::~Metric() {
	for (int i=0 ; i<espace.get_nbr_domains() ; i++) {
		if (p_met_cov[i]!=0x0)
			delete p_met_cov[i] ;
		p_met_cov[i] = 0x0 ;
		if (p_met_con[i]!=0x0)
			delete p_met_con[i] ;
		p_met_con[i] = 0x0 ;
		if (p_christo[i]!=0x0)
			delete p_christo[i] ;
		p_christo[i] = 0x0 ;		
		if (p_riemann[i]!=0x0)
			delete p_riemann[i] ;
		p_riemann[i] = 0x0 ;
		if (p_ricci_tensor[i]!=0x0)
			delete p_ricci_tensor[i] ;
		p_ricci_tensor[i] = 0x0 ;
		if (p_ricci_scalar[i]!=0x0)
			delete p_ricci_scalar[i] ;
		p_ricci_scalar[i] = 0x0 ;
		if (p_dirac[i]!=0x0)
			delete p_dirac[i] ;
		p_dirac[i] = 0x0 ;
		if (p_det_cov[i]!=0x0)
			delete p_det_cov[i] ;
		p_det_cov[i] = 0x0 ;
	}
	delete [] p_met_con ;
	delete [] p_met_cov ;	
	delete [] p_christo ;
	delete [] p_riemann ;
	delete [] p_ricci_tensor ;
	delete [] p_ricci_scalar ;
	delete [] p_dirac ;
	delete [] p_det_cov ;
	
}

void Metric::update(int i) {
  
	if (type_tensor == CON) {
	  if (p_met_con[i]!=0x0)
		  compute_con(i) ;
	  if (p_met_cov[i]!=0x0)
		compute_cov(i) ;
	}
	else {	  
	  if (p_met_cov[i]!=0x0)
		compute_cov(i) ;
	  if (p_met_con[i]!=0x0)
		  compute_con(i) ;
	}
	      if (p_det_cov[i]!=0x0)
			compute_det_cov(i) ;
		if (p_christo[i]!=0x0)
			compute_christo(i) ;		
		if (p_riemann[i]!=0x0)
			compute_riemann(i) ;
		if (p_dirac[i]!=0x0)
			compute_dirac(i) ;
		if (p_ricci_tensor[i]!=0x0)
			compute_ricci_tensor(i) ;
		if (p_ricci_scalar[i]!=0x0)
			compute_ricci_scalar(i) ;
			
}

void Metric::update() {
    for (int i=0 ; i<espace.get_nbr_domains() ; i++)
      update(i) ;
}

const Term_eq* Metric::give_term (int dd, int type) const {
	assert ((dd>=0) && (dd<espace.get_nbr_domains())) ;
	switch (type) {
		case COV :
			if (p_met_cov[dd]==0x0)
				compute_cov (dd) ;
			return p_met_cov[dd] ;
		case CON :
			if (p_met_con[dd]==0x0)
				compute_con (dd) ;
			return p_met_con[dd] ;
		default :
			cerr << "Unknown type of indice in Metric::give_term" << endl ;
			abort() ;
		}
}

int Metric::give_type (int dd) const {
  if (p_met_cov[dd]==0x0)
    compute_cov(dd) ;
  return (p_met_cov[dd]->val_t->get_basis().get_basis(dd)) ;
}

const Term_eq* Metric::give_christo (int dd) const {
	assert ((dd>=0) && (dd<espace.get_nbr_domains())) ;
	if (p_christo[dd]==0x0)
		compute_christo(dd) ;
	return p_christo[dd] ;
}

const Term_eq* Metric::give_riemann (int dd) const {
	assert ((dd>=0) && (dd<espace.get_nbr_domains())) ;
	if (p_riemann[dd]==0x0)
		compute_riemann(dd) ;
	return p_riemann[dd] ;
}

const Term_eq* Metric::give_ricci_tensor (int dd) const {
	assert ((dd>=0) && (dd<espace.get_nbr_domains())) ;
	if (p_ricci_tensor[dd]==0x0)
		compute_ricci_tensor(dd) ;
	return p_ricci_tensor[dd] ;
}

const Term_eq* Metric::give_ricci_scalar (int dd) const {
	assert ((dd>=0) && (dd<espace.get_nbr_domains())) ;
	if (p_ricci_scalar[dd]==0x0)
		compute_ricci_scalar(dd) ;
	return p_ricci_scalar[dd] ;
}

const Term_eq* Metric::give_dirac (int dd) const {
	assert ((dd>=0) && (dd<espace.get_nbr_domains())) ;
	if (p_dirac[dd]==0x0)
		compute_dirac(dd) ;
	return p_dirac[dd] ;
}

const Term_eq* Metric::give_det_cov (int dd) const {
	assert ((dd>=0) && (dd<espace.get_nbr_domains())) ;
	if (p_det_cov[dd]==0x0)
		compute_det_cov(dd) ;
	return p_det_cov[dd] ;
}


const Metric* Metric::get_background() const {
	cerr << "No background metric defined" << endl ; 
	abort() ;
}

Term_eq Metric::derive_partial (int type_der, char ind_der, const Term_eq& so) const {

	int dom = so.get_dom() ;

	if (type_tensor == CON) {
	  if (p_met_con[dom]==0x0)
		  compute_con(dom) ;
	  if (p_met_cov[dom]==0x0)
		compute_cov(dom) ;
	}
	else {	  
	  if (p_met_cov[dom]==0x0)
		compute_cov(dom) ;
	  if (p_met_con[dom]==0x0)
		  compute_con(dom) ;
	}
	
	// so must be tensor :
	if (so.get_type_data()!=TERM_T) {
		cerr << "Metric::derive partial only defined for tensor data" << endl ;
		abort() ;
	}

	// Check that the triad is good :
	if ((so.val_t->get_valence()>0) && (so.val_t->get_basis().get_basis(dom) != CARTESIAN_BASIS)) {
		cerr << "Metric::derive_partial only defined for tensor on cartesian triad" << endl ;
		abort() ;
	}
	
	
	bool doder = ((so.der_t==0x0) || (p_met_cov[dom]->der_t==0x0) || (p_met_cov[dom]->der_t==0x0)) ? false : true ;

	assert ((type_der==COV) || (type_der==CON)) ;
	int val_res = so.val_t->get_valence() + 1 ;

	bool doname = true ;
	if (so.val_t->get_valence()>0)
		if (!so.val_t->is_name_affected()) 
			doname = false;

	Array<int> type_ind (val_res) ;
	type_ind.set(0) = COV ;
	for (int i=1 ; i<val_res ; i++)
		type_ind.set(i) = so.val_t->get_index_type(i-1) ;

	// Need for summation ?
	bool need_sum = false ;
	if (doname)
		for (int i=1 ; i<val_res ; i++)
			if (ind_der== so.val_t->get_name_ind()[i-1])
				need_sum = true ;

	// Tensor for val
	Tensor auxi_val (espace, val_res, type_ind, so.val_t->get_basis()) ;
	// Set the names of the indices :
	if (doname) {
		auxi_val.set_name_affected() ;
		auxi_val.set_name_ind(0, ind_der) ;
		for (int i=1 ; i<val_res ; i++)
			auxi_val.set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;
	}

	//Loop on the components :
	Index pos_auxi(auxi_val) ;
	Index pos_so (*so.val_t) ;
	do {
		for (int i=0 ; i<val_res-1 ; i++)
			pos_so.set(i) = pos_auxi(i+1) ;
		auxi_val.set(pos_auxi).set_domain(dom) = (*so.val_t)(pos_so)(dom).der_abs(pos_auxi(0)+1) ;
	}
	while (pos_auxi.inc()) ;

	if (!doder) {
		// No need for derivative :
		Term_eq auxi (dom, auxi_val) ;
		// If derive contravariant : manipulate first indice :
		if (type_der==CON) 
			manipulate_ind (auxi, 0) ;
	
		if (!need_sum)
			return auxi ;
		else
			return Term_eq (dom, auxi.val_t->do_summation_one_dom(dom)) ;
	}
	else {
		// Need to compute the derivative :
		// Tensor for der
		Tensor auxi_der (espace, val_res, type_ind, so.der_t->get_basis()) ;
		// Set the names of the indices :
		auxi_der.set_name_affected() ;
		auxi_der.set_name_ind(0, ind_der) ;
		for (int i=1 ; i<val_res ; i++)
			auxi_der.set_name_ind(i, so.der_t->get_name_ind()[i-1]) ;

		//Loop on the components :
		Index pos_auxi_der(auxi_der) ;
		do {
			for (int i=0 ; i<val_res-1 ; i++)
				pos_so.set(i) = pos_auxi_der(i+1) ;
			auxi_der.set(pos_auxi_der).set_domain(dom) = (*so.der_t)(pos_so)(dom).der_abs(pos_auxi_der(0)+1) ;
		}
		while (pos_auxi_der.inc()) ;

		// Need for derivative :
		Term_eq auxi (dom, auxi_val, auxi_der) ;  
		// If derive contravariant : manipulate first indice :
		if (type_der==CON)
		    manipulate_ind (auxi, 0) ;
		
		if (!need_sum)
			return auxi ;
		else
		  return Term_eq (dom, auxi.val_t->do_summation_one_dom(dom), 
							auxi.der_t->do_summation_one_dom(dom)) ;
	}

}

Term_eq Metric::derive (int type_der, char ind_der, const Term_eq& so) const {

	// The partial derivative part
	Term_eq res (derive_partial (type_der, ind_der, so)) ;

	// Add the part containing the Christoffel :
	//Must find a name for summation on Christofel :
	bool found = false ;
	int start = 97 ;
	do {
		bool same = false ;
		if (ind_der==char(start))
			same = true ;
		for (int i=0 ; i<so.val_t->get_valence() ; i++)
			if (so.val_t->get_name_ind()[i]==char(start))
				same = true ;
		if (!same)
			found = true ;
		else 
			start ++ ;
	}
	while ((!found) && (start<123)) ;
	if (!found) {
		cerr << "Trouble with indices in derive (you are not using tensors of order > 24, are you ?)" << endl ;
		abort() ;
	}
	char name_sum = char(start) ;

	int dd = so.get_dom() ;
	if (p_christo[dd]==0x0)
		compute_christo(dd) ;
	bool doder = ((so.der_t==0x0) || (p_christo[dd]->der_t==0x0)) ? false : true ;

	//loop on the components of so :
	for (int cmp=0 ; cmp<so.val_t->get_valence() ; cmp ++) {

		int genre_indice = so.val_t->get_index_type(cmp) ;
		
		//Affecte names on the Christoffels :
		p_christo[dd]->val_t->set_name_affected() ;
		p_christo[dd]->val_t->set_name_ind(0, ind_der) ;
		if (genre_indice==COV) {
			p_christo[dd]->val_t->set_name_ind(1, so.val_t->get_name_ind()[cmp]) ;
			p_christo[dd]->val_t->set_name_ind(2, name_sum) ;
		}
		else {
			p_christo[dd]->val_t->set_name_ind(2, so.val_t->get_name_ind()[cmp]) ;
			p_christo[dd]->val_t->set_name_ind(1, name_sum) ;
		}

		if (doder) {
			p_christo[dd]->der_t->set_name_affected() ;
			p_christo[dd]->der_t->set_name_ind(0, ind_der) ;
			if (genre_indice==COV) {
				p_christo[dd]->der_t->set_name_ind(1, so.der_t->get_name_ind()[cmp]) ;
				p_christo[dd]->der_t->set_name_ind(2, name_sum) ;
			}
			else {
				p_christo[dd]->der_t->set_name_ind(2, so.der_t->get_name_ind()[cmp]) ;
				p_christo[dd]->der_t->set_name_ind(1, name_sum) ;
			}
		}

		Term_eq auxi_christ (*p_christo[dd]) ;

		if (type_der==CON) 
			manipulate_ind(auxi_christ, 0) ;

		// Check if one inner summation is needed :
		bool need_sum = false ;
		char* ind = p_christo[dd]->val_t->get_name_ind() ;
		if ((ind[0]==ind[2]) || (ind[1]==ind[2]) || (ind[0]==ind[1]))
			need_sum = true ;
				
		Term_eq* christ ;
		if (need_sum)
			if (!doder)
				christ = new Term_eq (dd, auxi_christ.val_t->do_summation_one_dom(dd)) ;
			else
				christ = new Term_eq (dd, auxi_christ.val_t->do_summation_one_dom(dd), 
					auxi_christ.der_t->do_summation_one_dom(dd)) ;
		else 
			christ = new Term_eq (auxi_christ) ;

		// Affecte names on the field :
		Term_eq copie (so) ;
		copie.val_t->set_name_affected() ;
		for (int i=0 ; i<so.val_t->get_valence() ; i++)
			if (i!=cmp)
				copie.val_t->set_name_ind(i, so.val_t->get_name_ind()[i]) ;
			else
				copie.val_t->set_name_ind(i, name_sum) ;
		if (doder) {
			copie.der_t->set_name_affected() ;
			for (int i=0 ; i<so.der_t->get_valence() ; i++)
				if (i!=cmp)
					copie.der_t->set_name_ind(i, so.der_t->get_name_ind()[i]) ;
				else
					copie.der_t->set_name_ind(i, name_sum) ;
		}

		Term_eq part_christo ((*christ)*copie) ;
		delete christ ;

		if (genre_indice==CON)
			res = res + part_christo ;
		else
			res = res - part_christo ;
	}
	return res ;
}

Term_eq Metric::derive_flat (int, char, const Term_eq&) const {
  cerr << "derive_flat not implemented for this type of metric (yet)" << endl ;
  abort() ;
}

void Metric::compute_christo (int dd) const {

	if (type_tensor == CON) {
	  if (p_met_con[dd]==0x0)
		  compute_con(dd) ;
	  if (p_met_cov[dd]==0x0)
		compute_cov(dd) ;
	}
	else {	  
	  if (p_met_cov[dd]==0x0)
		compute_cov(dd) ;
	  if (p_met_con[dd]==0x0)
		  compute_con(dd) ;
	}

	// Get the conformal factor
	bool doder = ((p_met_con[dd]->der_t==0x0) || (p_met_cov[dd]->der_t==0x0)) ? false : true ;

	Array<int> type_ind (3) ;
	type_ind.set(0) = COV ; type_ind.set(1) = COV ; type_ind.set(2) = CON ;
	Tensor res_val (espace, 3, type_ind, p_met_con[dd]->val_t->get_basis()) ;
	Tensor res_der (espace, 3, type_ind, p_met_con[dd]->val_t->get_basis()) ;
   
	res_val = 0 ;
	res_der = 0 ;
	
	Index pos (res_val) ;
	do {
		Val_domain cmpval(espace.get_domain(dd)) ;
		cmpval = 0 ;
		
		for (int l=1 ; l<=espace.get_ndim() ; l++)
			cmpval += 0.5*(
			(*p_met_con[dd]->val_t)(pos(2)+1,l)(dd)*((*p_met_cov[dd]->val_t)(pos(1)+1,l)(dd).der_abs(pos(0)+1) + 
				(*p_met_cov[dd]->val_t)(pos(0)+1,l)(dd).der_abs(pos(1)+1) - (*p_met_cov[dd]->val_t)(pos(1)+1,pos(0)+1)(dd).der_abs(l))) ;

		res_val.set(pos).set_domain(dd) = cmpval ;	

		if (doder) {
		    Val_domain cmpder(espace.get_domain(dd)) ;
		    cmpder = 0 ;
		    
		    for (int l=1 ; l<=espace.get_ndim() ; l++)
			cmpder +=  0.5*(
			(*p_met_con[dd]->der_t)(pos(2)+1,l)(dd)*((*p_met_cov[dd]->val_t)(pos(1)+1,l)(dd).der_abs(pos(0)+1) + 
				(*p_met_cov[dd]->val_t)(pos(0)+1,l)(dd).der_abs(pos(1)+1) - (*p_met_cov[dd]->val_t)(pos(1)+1,pos(0)+1)(dd).der_abs(l))) 
			+ 0.5*(
			(*p_met_con[dd]->val_t)(pos(2)+1,l)(dd)*((*p_met_cov[dd]->der_t)(pos(1)+1,l)(dd).der_abs(pos(0)+1) + 
				(*p_met_cov[dd]->der_t)(pos(0)+1,l)(dd).der_abs(pos(1)+1) - (*p_met_cov[dd]->der_t)(pos(1)+1,pos(0)+1)(dd).der_abs(l))) ;
		
		    res_der.set(pos).set_domain(dd) = cmpder ;
		}
	}
	while (pos.inc()) ;
	
	if (!doder) {
		if (p_christo[dd]==0x0)
			p_christo[dd] = new Term_eq(dd, res_val) ;
		else
			*p_christo[dd] = Term_eq (dd, res_val) ;
	}
	else {
		if (p_christo[dd]==0x0)
			p_christo[dd] = new Term_eq(dd, res_val, res_der) ;
		else
			*p_christo[dd] = Term_eq(dd, res_val, res_der) ;
	}
}

void Metric::compute_riemann (int dd) const {
	// Need christoffels
	if (p_christo[dd]==0x0)
		compute_christo(dd) ;

	// Get the conformal factor
	bool doder = ((p_met_con[dd]->der_t==0x0) || (p_met_cov[dd]->der_t==0x0)) ? false : true ;
	Array<int> indices (4) ;
	indices.set(0) = CON ; indices.set(1) = COV ; indices.set(2) = COV ; indices.set(3) = COV ;
	Tensor res_val (espace, 4, indices, p_met_con[dd]->val_t->get_basis()) ;
	Tensor res_der (espace, 4, indices, p_met_con[dd]->val_t->get_basis()) ;

	Index pos (res_val) ;
	do {
		Val_domain cmpval (espace.get_domain(dd)) ;

		cmpval = (*p_christo[dd]->val_t)(pos(1)+1,pos(3)+1,pos(0)+1)(dd).der_abs(pos(2)+1) - 
			(*p_christo[dd]->val_t)(pos(1)+1, pos(2)+1, pos(0)+1)(dd).der_abs(pos(3)+1) ;

		for (int m=1 ; m<=espace.get_ndim() ; m++)
			cmpval += (*p_christo[dd]->val_t)(pos(2)+1,m, pos(0)+1)(dd)*(*p_christo[dd]->val_t)(pos(1)+1,pos(3)+1,m)(dd) 
						- (*p_christo[dd]->val_t)(pos(3)+1,m,pos(0)+1)(dd)*(*p_christo[dd]->val_t)(pos(1)+1,pos(2)+1,m)(dd) ;
		res_val.set(pos).set_domain(dd) = cmpval ;
		
		if (doder) {
		  Val_domain cmpder(espace.get_domain(dd)) ;
		  
		  cmpder = (*p_christo[dd]->der_t)(pos(1)+1,pos(3)+1,pos(0)+1)(dd).der_abs(pos(2)+1) - 
			(*p_christo[dd]->der_t)(pos(1)+1, pos(2)+1, pos(0)+1)(dd).der_abs(pos(3)+1) ;
		
		for (int m=1 ; m<=espace.get_ndim() ; m++)
			cmpder += (*p_christo[dd]->der_t)(pos(2)+1,m, pos(0)+1)(dd)*(*p_christo[dd]->val_t)(pos(1)+1,pos(3)+1,m)(dd) 
			+   (*p_christo[dd]->val_t)(pos(2)+1,m, pos(0)+1)(dd)*(*p_christo[dd]->der_t)(pos(1)+1,pos(3)+1,m)(dd) 
			- (*p_christo[dd]->der_t)(pos(3)+1,m,pos(0)+1)(dd)*(*p_christo[dd]->val_t)(pos(1)+1,pos(2)+1,m)(dd) 
			- (*p_christo[dd]->val_t)(pos(3)+1,m,pos(0)+1)(dd)*(*p_christo[dd]->der_t)(pos(1)+1,pos(2)+1,m)(dd) ;

		res_der.set(pos).set_domain(dd) = cmpder ;
		}
	}
	while (pos.inc()) ;

	if (!doder) {
		if (p_riemann[dd]==0x0)
			p_riemann[dd] = new Term_eq(dd, res_val) ;
		else
			*p_riemann[dd] = Term_eq (dd, res_val) ;
	}
	else {
		if (p_riemann[dd]==0x0)
			p_riemann[dd] = new Term_eq(dd, res_val, res_der) ;
		else
			*p_riemann[dd] = Term_eq(dd, res_val, res_der) ;
	}	
}


void Metric::compute_ricci_tensor (int dd) const {
	// Need christoffels
	if (p_christo[dd]==0x0)
		compute_christo(dd) ;

	// Get the conformal factor
	bool doder = ((p_met_con[dd]->der_t==0x0) || (p_met_cov[dd]->der_t==0x0)) ? false : true ;
	Tensor res_val (espace, 2, COV, p_met_con[dd]->val_t->get_basis()) ;
	Tensor res_der (espace, 2, COV, p_met_con[dd]->val_t->get_basis()) ;

	Index pos (res_val) ;
	do {
		Val_domain cmpval (espace.get_domain(dd)) ;
		cmpval = 0 ;
		
		for (int k=1 ; k<=espace.get_ndim() ; k++) {
			cmpval += (*p_christo[dd]->val_t)(pos(0)+1,pos(1)+1,k)(dd).der_abs(k) - (*p_christo[dd]->val_t)(pos(1)+1, k, k)(dd).der_abs(pos(0)+1) ;
			for (int l=1 ; l<=espace.get_ndim() ; l++)
				cmpval += (*p_christo[dd]->val_t)(k,l,l)(dd)*(*p_christo[dd]->val_t)(pos(0)+1,pos(1)+1,k)(dd) 
						- (*p_christo[dd]->val_t)(pos(0)+1,k,l)(dd)*(*p_christo[dd]->val_t)(pos(1)+1,l,k)(dd) ;
		}
		
		res_val.set(pos).set_domain(dd) = cmpval ;	

		if (doder) {
		    Val_domain cmpder (espace.get_domain(dd)) ;
		    cmpder = 0 ;
		    
		for (int k=1 ; k<=espace.get_ndim() ; k++) {
			cmpder += (*p_christo[dd]->der_t)(pos(0)+1,pos(1)+1,k)(dd).der_abs(k) - (*p_christo[dd]->der_t)(pos(1)+1, k, k)(dd).der_abs(pos(0)+1) ;
			for (int l=1 ; l<=espace.get_ndim() ; l++)
				cmpder += (*p_christo[dd]->der_t)(k,l,l)(dd)*(*p_christo[dd]->val_t)(pos(0)+1,pos(1)+1,k)(dd) 
					  +  (*p_christo[dd]->val_t)(k,l,l)(dd)*(*p_christo[dd]->der_t)(pos(0)+1,pos(1)+1,k)(dd) 
					  - (*p_christo[dd]->der_t)(pos(0)+1,k,l)(dd)*(*p_christo[dd]->val_t)(pos(1)+1,l,k)(dd)
					  - (*p_christo[dd]->val_t)(pos(0)+1,k,l)(dd)*(*p_christo[dd]->der_t)(pos(1)+1,l,k)(dd) ;
		}
		res_der.set(pos).set_domain(dd) = cmpder ;
		}
	}
	while (pos.inc()) ;

	if (!doder) {
		if (p_ricci_tensor[dd]==0x0)
			p_ricci_tensor[dd] = new Term_eq(dd, res_val) ;
		else
			*p_ricci_tensor[dd] = Term_eq (dd, res_val) ;
	}
	else {
		if (p_ricci_tensor[dd]==0x0)
			p_ricci_tensor[dd] = new Term_eq(dd, res_val, res_der) ;
		else
			*p_ricci_tensor[dd] = Term_eq(dd, res_val, res_der) ;
	}

      
	
}

void Metric::compute_ricci_scalar (int dd) const {

	// Need that
	if (type_tensor == CON) {
	  if (p_met_con[dd]==0x0)
		  compute_con(dd) ;
	  if (p_met_cov[dd]==0x0)
		compute_cov(dd) ;
	}
	else {	  
	  if (p_met_cov[dd]==0x0)
		compute_cov(dd) ;
	  if (p_met_con[dd]==0x0)
		  compute_con(dd) ;
	}
	if (p_ricci_tensor[dd]==0x0)
		compute_ricci_tensor(dd) ;

	bool doder = ((p_met_con[dd]->der_t==0x0) || (p_ricci_tensor[dd]->der_t==0x0)) ? false : true ;
	Scalar res_val (espace) ;
	Scalar res_der (espace) ;

	Val_domain cmpval (espace.get_domain(dd)) ;
	cmpval = 0 ;
	
	for (int i=1 ; i<=espace.get_ndim() ; i++)
		for (int j=1 ; j<=espace.get_ndim() ; j++)
			cmpval += (*p_met_con[dd]->val_t)(i,j)(dd)*(*p_ricci_tensor[dd]->val_t)(i,j)(dd) ;
	res_val.set_domain(dd) = cmpval ;
	
	if (doder) {
	  Val_domain cmpder (espace.get_domain(dd)) ;
	  cmpder = 0 ;
	  
	  for (int i=1 ; i<=espace.get_ndim() ; i++)
		for (int j=1 ; j<=espace.get_ndim() ; j++)
			cmpder += (*p_met_con[dd]->val_t)(i,j)(dd)*(*p_ricci_tensor[dd]->der_t)(i,j)(dd) 
				  + (*p_met_con[dd]->der_t)(i,j)(dd)*(*p_ricci_tensor[dd]->val_t)(i,j)(dd) ;
	  res_der.set_domain(dd) = cmpder ;
	}
	
	if (!doder) {
		if (p_ricci_scalar[dd]==0x0)
			p_ricci_scalar[dd] = new Term_eq(dd, res_val) ;
		else
			*p_ricci_scalar[dd] = Term_eq (dd, res_val) ;
	}
	else {
		if (p_ricci_scalar[dd]==0x0)
			p_ricci_scalar[dd] = new Term_eq(dd, res_val, res_der) ;
		else
			*p_ricci_scalar[dd] = Term_eq(dd, res_val, res_der) ;
	}
}

void Metric::manipulate_ind(Term_eq& so, int ind) const {
	
	int dd = so.get_dom() ;
	int valence = so.get_val_t().get_valence() ;
	int type_start = so.get_val_t().get_index_type (ind) ;
	
	if (type_tensor == CON) {
	  if (p_met_con[dd]==0x0)
		  compute_con(dd) ;
	  if (p_met_cov[dd]==0x0)
		compute_cov(dd) ;
	}
	else {	  
	  if (p_met_cov[dd]==0x0)
		compute_cov(dd) ;
	  if (p_met_con[dd]==0x0)
		  compute_con(dd) ;
	}
	
	bool doder = false ;
	if ((type_start==COV) && (p_met_con[dd]->der_t!=0x0) && (so.der_t!=0x0))
		doder = true ;
	if ((type_start==CON) && (p_met_cov[dd]->der_t!=0x0) && (so.der_t!=0x0))
		doder = true ;

	Array<int> type_res (valence) ;
	for (int i=0 ; i<valence ; i++)
		type_res.set(i) = (i==ind) ? - so.get_val_t().get_index_type(i) : so.get_val_t().get_index_type(i) ;
	
	Tensor val_res (espace, valence, type_res, p_met_con[dd]->val_t->get_basis()) ;
	Tensor val_der (espace, valence, type_res, p_met_con[dd]->val_t->get_basis()) ;
	
	Index pos (val_res) ;
	do {
		Val_domain cmpval (espace.get_domain(dd)) ;
		cmpval = 0 ;
		
		for (int k=0 ; k<espace.get_ndim() ; k++) {
			Index copie(pos) ;
			copie.set(ind) = k ;
			if (type_start==COV)
				cmpval += (*p_met_con[dd]->val_t)(pos(ind)+1, k+1)(dd) * (*so.val_t)(copie)(dd);
			else
				cmpval +=  (*p_met_cov[dd]->val_t)(pos(ind)+1, k+1)(dd) * (*so.val_t)(copie)(dd) ;
		}
		
		val_res.set(pos).set_domain(dd) = cmpval ;
		
		if (doder) {
		  Val_domain cmpder(espace.get_domain(dd)) ;
		  cmpder = 0 ;
		  for (int k=0 ; k<espace.get_ndim() ; k++) {
			Index copie(pos) ;
			copie.set(ind) = k ;
			if (type_start==COV)
				cmpder += (*p_met_con[dd]->der_t)(pos(ind)+1, k+1)(dd) * (*so.val_t)(copie)(dd) 
					+(*p_met_con[dd]->val_t)(pos(ind)+1, k+1)(dd) * (*so.der_t)(copie)(dd);
			else
				cmpder +=  (*p_met_cov[dd]->der_t)(pos(ind)+1, k+1)(dd) * (*so.val_t)(copie)(dd) 
					+  (*p_met_cov[dd]->val_t)(pos(ind)+1, k+1)(dd) * (*so.der_t)(copie)(dd) ;	
				
		}
		val_der.set(pos).set_domain(dd) = cmpder ;
		}
	}
	while (pos.inc()) ;

	// Put the names :
	if (so.get_val_t().is_name_affected()) {
			val_res.set_name_affected() ;
			for (int i=0 ; i<valence ; i++)
				val_res.set_name_ind(i, so.val_t->get_name_ind()[i]) ;
	}

	if ((doder) && (so.get_der_t().is_name_affected())) {
			val_der.set_name_affected() ;
			for (int i=0 ; i<valence ; i++)
				val_der.set_name_ind(i, so.der_t->get_name_ind()[i]) ;
	}

	so.set_val_t()->set_index_type (ind) *= -1 ;
	if (so.set_der_t() !=0x0)
		so.set_der_t()->set_index_type (ind) *= -1 ;
	
	delete so.val_t ;
	so.val_t = new Tensor (val_res) ;
	delete so.der_t ;
	so.der_t = (doder) ? new Tensor (val_der) : 0x0 ;
}

void Metric::compute_cov (int) const {
	cerr << "compute_cov not implemented for this type of metric" << endl ;
	abort() ;
}

void Metric::compute_con (int) const {
	cerr << "compute_con not implemented for this type of metric" << endl ;
	abort() ;
}

void Metric::compute_dirac (int) const {
    cerr << "Compute Dirac not implemented for this type of metric" << endl ;
    abort() ;
}

void Metric::compute_det_cov (int) const {
    cerr << "Compute determinant of covariant not implemented for this type of metric" << endl ;
    abort() ;
}}
