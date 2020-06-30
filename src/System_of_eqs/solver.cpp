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

#include "system_of_eqs.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "metric.hpp"
#include "exceptions.hpp"
namespace Kadath {
Array<double> System_of_eqs::check_equations() {

	
	Array<double> sec (sec_member()) ;
	Array<double> errors (neq_int + neq) ;
	errors = 0 ;
      
	int pos = 0 ;
	for (int i=0 ; i<neq_int ; i++) {
	  errors.set(i) = fabs(sec(pos)) ;
	  pos ++ ;
	}

	for (int i=0 ; i<neq ; i++) {
	    double max = 0 ;
	    for (int j=0 ; j<eq[i]->get_n_cond_tot() ; j++) {
	      if (fabs(sec(pos)) > max)
		max = fabs(sec(pos)) ;
	      pos ++ ;
	    }
	    errors.set(neq_int+i) = max ;
	}
      return errors ;
}

void System_of_eqs::compute_nbr_of_conditions()
{

    if (met!=0x0)
        for (int d=dom_min ; d<=dom_max ; d++)
            met->update(d) ;
    for (int i=0 ; i<ndef ; i++)
        def[i]->compute_res() ;

    int conte = 0 ;
    for (int i=0 ; i<neq ; i++)
        eq[i]->apply(conte, results.set_data()) ;

    // Need to assert the size :
    if (nbr_conditions==-1) {
        nbr_conditions = 0 ;
        for (int i=0 ; i<neq_int ; i++)
            nbr_conditions ++ ;
        for (int i=0 ; i<neq ; i++)
            nbr_conditions += eq[i]->get_n_cond_tot() ;
    }

}

Array<double> System_of_eqs::sec_member() {

    this->vars_to_terms() ;
    this->compute_nbr_of_conditions();
	// Computation of the second member itself :
	Array<double> res (nbr_conditions) ;
	res = 0 ;
	int conte = 0 ;
	int pos_res = 0 ;
	for (int i=0 ; i<neq_int ; i++) {
		res.set(pos_res) = eq_int[i]->get_val() ;
		pos_res ++ ;
	}

	for (int i=0 ; i<neq ; i++)
	    eq[i]->export_val(conte, results.set_data(), res, pos_res) ;
	return res ;
}

Array<double> System_of_eqs::do_JX (const Array<double>& xx) {

	xx_to_ders(xx) ;	
	if (met!=0x0)  
	  for (int d=dom_min ; d<=dom_max ; d++)
		met->update(d) ;
	// Delete the definitions 
	for (int i=0 ; i<ndef ; i++)
		def[i]->compute_res() ;

	int conte = 0 ;
	for (int i=0 ; i<neq ; i++)  
	    eq[i]->apply(conte, results.set_data()) ;
	
	if (nbr_conditions==-1) {
		cerr << "Number of conditions unknown ; call sec_member first" << endl ;
		abort() ;
	}
	
	Array<double> res (nbr_conditions) ;
	conte = 0 ;
	int pos_res = 0 ;
	res = 0 ;
	for (int i=0 ; i<neq_int ; i++) {
		res.set(pos_res) = eq_int[i]->get_der() ;
		pos_res ++ ;
	}

	for (int i=0 ; i<neq ; i++)
	    eq[i]->export_der(conte, results.set_data(), res, pos_res) ;
	return res ;
}

Array<double> System_of_eqs::do_col_J (int cc) {

	assert ((cc>=0) && (cc<nbr_unknowns)) ;

	// Affecte nterms derivatives :
	int conte = 0 ;
	int zedom = -1 ;
	bool is_var_double = false ;
	Array<int> zedoms (2) ;
	zedoms = - 1 ;
	// Variable Domains :
	espace.affecte_coef_to_variable_domains(conte, cc, zedoms) ;
	if (zedoms(0)!=-1)
	  update_terms_from_variable_domains(zedoms) ;
	else {
	  for (int i=0 ; i<nterm_cst ; i++)
	    cst[i]->set_der_zero() ;
	  for (int i=0 ; i<nterm ; i++)
	    term[i]->set_der_zero() ;
	}
	
	// Double 
	for (int i=0 ; i<nvar_double ; i++) {
		if (conte==cc) {
			for (int dd=dom_min ; dd<=dom_max ; dd++)
				term_double[i*ndom+(dd-dom_min)]->set_der_d(1.) ;
			is_var_double = true ;
		}
		else
			for (int dd=dom_min ; dd<=dom_max ; dd++)
				term_double[i*ndom+(dd-dom_min)]->set_der_d(0.) ;
		conte ++ ;
	}

	// Fields
	for (int i=0 ; i<nterm ; i++) {
		int dom = term[i]->get_dom() ;
		Tensor auxi (term[i]->get_val_t(), false) ;
		try {
            espace.get_domain(dom)->affecte_tau_one_coef(auxi, dom, cc, conte);
        }
		catch(Unknown_base_error & e) {
		    std::cerr << "Error in System_of_eqs[mpi_proc_rank=" << mpi_proc_rank << "]::do_col_J(cc = " << cc
                     << ")\n";
		    std::cerr << "affecte_tau_one_coef raised the following exception : " << e.what() << "\n";
		    std::cerr << "while calling espace.get_domain(dom)->affecte_tau_one_coef(auxi, dom, cc, conte)\n";
		    std::cerr << "with auxi=term[i]->get_val_t(), where i="<< i << " and *term[i] = \n" << *term[i];
		    std::cerr << "\n dom=" << dom << "and conte=" << conte << std::endl;
		    abort();
		}
		for (int j=0 ; j<auxi.get_n_comp() ; j++) {
			// Si la base n'est pas affectee on la met
			if ((!auxi(auxi.indices(j))(dom).check_if_zero()) && (!auxi(auxi.indices(j))(dom).get_base().is_def()))
				auxi.set(auxi.indices(j)).set_domain(dom).set_base() = term[i]->get_val_t()(auxi.indices(j))(dom).get_base() ;
			if ((zedom==-1) && (!auxi(auxi.indices(j))(dom).check_if_zero()))
				zedom = dom ;
		}

		term[i]->set_der_t(auxi + term[i]->get_der_t()) ;
	}
	
	// Delete the metric derivative terms :
	if (met!=0x0)  
	  for (int d=dom_min ; d<=dom_max ; d++)
		met->update(d) ;
	  
	// Delete the definitions 
	for (int i=0 ; i<ndef ; i++)
		def[i]->compute_res() ;
	
	conte = 0 ;
	for (int i=0 ; i<neq ; i++) {
		if ((is_var_double) || (eq[i]->take_into_account(zedom)) || (eq[i]->take_into_account(zedoms(0))) || (eq[i]->take_into_account(zedoms(1)))) {
			eq[i]->apply(conte, results.set_data()) ;//NOT THREAD SAFE (but both eq and results are own by System_of_eqs)
		}
		else
			conte += eq[i]->n_ope ;
	}
		
	if (nbr_conditions==-1) {
		cerr << "Number of conditions unknown ; call sec_member first" << endl ;
		abort() ;
	}
	
	Array<double> res (nbr_conditions) ;
	conte = 0 ;
	int pos_res = 0 ;
	res = 0 ;
	for (int i=0 ; i<neq_int ; i++) {
		res.set(pos_res) = eq_int[i]->get_der() ;
		pos_res ++ ;
	} 
	
	for (int i=0 ; i<neq ; i++) {
		if ((is_var_double) || (eq[i]->take_into_account(zedom))|| (eq[i]->take_into_account(zedoms(0))) || (eq[i]->take_into_account(zedoms(1)))) {
		  eq[i]->export_der(conte, results.set_data(), res, pos_res) ;
		}
		else
			{
			pos_res += eq[i]->n_cond_tot ;
			conte += eq[i]->n_ope ;
			}  
	}
	
	return res ;
}

void System_of_eqs::update_terms_from_variable_domains(const Array<int>& zedoms) {
  
	for (int i=0 ; i<nterm ; i++) {
	      int dom = term[i]->get_dom() ;
	      bool  todo = false ;
	      for (int d=0 ; d<zedoms.get_size(0) ; d++)
		  if (zedoms(d)==dom)
		      todo = true ;
	      if (todo)
		espace.get_domain(dom)->update_term_eq (term[i]) ;
	      else
		term[i]->set_der_zero() ;
	}
	  
	  for (int i=0 ; i<nterm_cst ; i++) 
	      if (cst[i]->get_type_data() == TERM_T) {
	      int dom = cst[i]->get_dom() ;
	      bool  todo = false ;
	      for (int d=0 ; d<zedoms.get_size(0) ; d++)
		  if (zedoms(d)==dom)
		      todo = true ;
	      if (todo)
		espace.get_domain(dom)->update_term_eq (cst[i]) ;
	      else
		cst[i]->set_der_zero() ;
	  } 
}}
