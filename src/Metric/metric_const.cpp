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
#include "metric_tensor.hpp"
#include "name_tools.hpp"
namespace Kadath {
Metric_const::Metric_const (Metric_tensor& met) : 
		Metric_general(met) {
}

Metric_const::Metric_const (const Metric_const& so) :
		Metric_general (*so.p_met) {
}

Metric_const::~Metric_const() {
}

void Metric_const::set_system (System_of_eqs& ss, const char* name_met) {

	syst = &ss ;

	if (ss.met!=0x0) {
		cerr << "Metric already set for the system" << endl ;
		abort() ;
	}

	// Position in the system :
        place_syst = ss.ndom*ss.ncst ;

	ss.add_cst (0x0, *p_met) ;

	ss.met = this ;
	ss.name_met = new char[LMAX] ;
	trim_spaces (ss.name_met, name_met) ;
}



void Metric_const::compute_cov (int dd) const {

	int dim = espace.get_ndim() ;
	if (dim!=3) {
		cerr << "Function only implemented for dimension 3" << endl ;
		abort() ;
	}

	int place = place_syst + (dd-syst->dom_min) ;

	// Right storage : simple copy.
	if (type_tensor==COV) {

		if (p_met_cov[dd]==0x0)
			p_met_cov[dd] = new Term_eq(*syst->cst[place]) ;
		else
			*p_met_cov[dd] = Term_eq(*syst->cst[place]) ;
	}
	else {
		Term_eq** res = new Term_eq* [p_met->get_n_comp()] ;
		
		
		Scalar val (espace) ;
		Val_domain cmpval(espace.get_domain(dd)) ;
	
		Val_domain detval (espace.get_domain(dd)) ;
		detval = (*syst->cst[place]->val_t)(1,1)(dd)*(*syst->cst[place]->val_t)(2,2)(dd)*(*syst->cst[place]->val_t)(3,3)(dd)
		      + (*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(2,3)(dd)*(*syst->cst[place]->val_t)(1,3)(dd)
		      + (*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(2,3)(dd)
		      - (*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(2,2)(dd)
		      - (*syst->cst[place]->val_t)(2,3)(dd)*(*syst->cst[place]->val_t)(2,3)(dd)*(*syst->cst[place]->val_t)(1,1)(dd)
		      - (*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(3,3)(dd)  ;

		// Compo 1 1
		cmpval = (*syst->cst[place]->val_t)(2,2)(dd)*(*syst->cst[place]->val_t)(3,3)(dd)
		      -(*syst->cst[place]->val_t)(2,3)(dd)*(*syst->cst[place]->val_t)(2,3)(dd) ;
		val.set_domain(dd) = cmpval / detval ;
		res[0] = new Term_eq (dd, val) ;
		
		  
		// Compo 1 2
		cmpval = (*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(2,3)(dd)
		      -(*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(3,3)(dd) ;
		val.set_domain(dd) = cmpval / detval ;
		res[1] = new Term_eq (dd, val) ;
		 
		// Compo 1 3
		cmpval = (*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(2,3)(dd)
		      -(*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(2,2)(dd) ;
		val.set_domain(dd) = cmpval / detval ;
		res[2] = new Term_eq (dd, val) ;
		
		
		// Compo 2 2
		cmpval = (*syst->cst[place]->val_t)(1,1)(dd)*(*syst->cst[place]->val_t)(3,3)(dd)
		      -(*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(1,3)(dd) ;
		val.set_domain(dd) = cmpval /detval ;
		res[3] = new Term_eq (dd, val) ;
		
		// Compo 2 3
		cmpval = (*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(1,2)(dd)
		      -(*syst->cst[place]->val_t)(1,1)(dd)*(*syst->cst[place]->val_t)(2,3)(dd) ;
		val.set_domain(dd) = cmpval / detval ;
		res[4] = new Term_eq (dd, val) ;
		
		
	      // Compo 3 3
		cmpval = (*syst->cst[place]->val_t)(1,1)(dd)*(*syst->cst[place]->val_t)(2,2)(dd)
		      -(*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(1,2)(dd) ;
		val.set_domain(dd) = cmpval /detval ;
		res[5] = new Term_eq (dd, val) ;
		
		
		// Value field :
		Metric_tensor resval (espace, COV, basis) ;
		resval.set(1,1) = res[0]->get_val_t() ;
		resval.set(1,2) = res[1]->get_val_t() ;
		resval.set(1,3) = res[2]->get_val_t() ;
		resval.set(2,2) = res[3]->get_val_t() ;
		resval.set(2,3) = res[4]->get_val_t() ;
		resval.set(3,3) = res[5]->get_val_t() ;

		Metric_tensor zero (espace, COV, basis) ;
		for (int i=1 ; i<=3 ; i++)
		  for (int j=i ; j<=3 ; j++)
		    zero.set(i,j) = 0 ;
		
		if (p_met_cov[dd]==0x0)
		    p_met_cov[dd] = new Term_eq(dd, resval, zero) ;
		else
		    *p_met_cov[dd] = Term_eq(dd, resval, zero) ;
		for (int i=0 ; i<p_met->get_n_comp() ; i++)
			delete res [i] ;
		delete [] res ;
	}
}


void Metric_const::compute_con (int dd) const {
	int dim = espace.get_ndim() ;
	if (dim!=3) {
		cerr << "Function only implemented for dimension 3" << endl ;
		abort() ;
	}

	int place = place_syst + (dd-syst->dom_min) ;

	// Right storage : simple copy.
	if (type_tensor==CON) {

		if (p_met_con[dd]==0x0)
			p_met_con[dd] = new Term_eq(*syst->cst[place]) ;
		else
			*p_met_con[dd] = Term_eq(*syst->cst[place]) ;
	}
	else {
	  
		Term_eq** res = new Term_eq* [p_met->get_n_comp()] ;
		
		
		Scalar val (espace) ;
		Val_domain cmpval(espace.get_domain(dd)) ;
		Val_domain detval (espace.get_domain(dd)) ;
		detval = (*syst->cst[place]->val_t)(1,1)(dd)*(*syst->cst[place]->val_t)(2,2)(dd)*(*syst->cst[place]->val_t)(3,3)(dd)
		      + (*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(2,3)(dd)*(*syst->cst[place]->val_t)(1,3)(dd)
		      + (*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(2,3)(dd)
		      - (*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(2,2)(dd)
		      - (*syst->cst[place]->val_t)(2,3)(dd)*(*syst->cst[place]->val_t)(2,3)(dd)*(*syst->cst[place]->val_t)(1,1)(dd)
		      - (*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(3,3)(dd)  ;

		// Compo 1 1
		cmpval = (*syst->cst[place]->val_t)(2,2)(dd)*(*syst->cst[place]->val_t)(3,3)(dd)
		      -(*syst->cst[place]->val_t)(2,3)(dd)*(*syst->cst[place]->val_t)(2,3)(dd) ;
		val.set_domain(dd) = cmpval /detval ;
		res[0] = new Term_eq (dd, val) ;
		  
		// Compo 1 2
		cmpval = (*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(2,3)(dd)
		      -(*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(3,3)(dd) ;
		val.set_domain(dd) = cmpval /detval;
		res[1] = new Term_eq (dd, val) ;
		  
		// Compo 1 3
		cmpval = (*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(2,3)(dd)
		      -(*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(2,2)(dd) ;
		val.set_domain(dd) = cmpval /detval;
		res[2] = new Term_eq (dd, val) ;
		
		// Compo 2 2
		cmpval = (*syst->cst[place]->val_t)(1,1)(dd)*(*syst->cst[place]->val_t)(3,3)(dd)
		      -(*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(1,3)(dd) ;
		val.set_domain(dd) = cmpval/detval ;
		res[3] = new Term_eq (dd, val) ;
		
		// Compo 2 3
		cmpval = (*syst->cst[place]->val_t)(1,3)(dd)*(*syst->cst[place]->val_t)(1,2)(dd)
		      -(*syst->cst[place]->val_t)(1,1)(dd)*(*syst->cst[place]->val_t)(2,3)(dd) ;
		val.set_domain(dd) = cmpval /detval;
		res[4] = new Term_eq (dd, val) ;
		
	      // Compo 3 3
		cmpval = (*syst->cst[place]->val_t)(1,1)(dd)*(*syst->cst[place]->val_t)(2,2)(dd)
		      -(*syst->cst[place]->val_t)(1,2)(dd)*(*syst->cst[place]->val_t)(1,2)(dd) ;
		val.set_domain(dd) = cmpval /detval;
		res[5] = new Term_eq (dd, val) ;
		
		// Value field :
		Metric_tensor resval (espace, CON, basis) ;
		resval.set(1,1) = res[0]->get_val_t() ;
		resval.set(1,2) = res[1]->get_val_t() ;
		resval.set(1,3) = res[2]->get_val_t() ;
		resval.set(2,2) = res[3]->get_val_t() ;
		resval.set(2,3) = res[4]->get_val_t() ;
		resval.set(3,3) = res[5]->get_val_t() ;
		
		Metric_tensor zero (espace, CON, basis) ;
		for (int i=1 ; i<=3 ; i++)
		  for (int j=i ; j<=3 ; j++)
		    zero.set(i,j) = 0 ;
		
		if (p_met_con[dd]==0x0)
			p_met_con[dd] = new Term_eq(dd, resval, zero) ;
		else
			*p_met_con[dd] = Term_eq(dd, resval, zero) ;
		for (int i=0 ; i<p_met->get_n_comp() ; i++)
			delete res [i] ;
		delete [] res ;
	}
}}
