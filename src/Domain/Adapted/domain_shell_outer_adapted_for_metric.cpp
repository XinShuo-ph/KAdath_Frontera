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
#include "utilities.hpp"
#include "adapted.hpp"
#include "point.hpp"
#include "array_math.cpp"
#include "term_eq.hpp"
#include "metric.hpp"
#include "scalar.hpp"

namespace Kadath {
Term_eq Domain_shell_outer_adapted::derive_flat_spher (int type_der, char ind_der, const Term_eq& so, const Metric* manipulator) const {
	assert ((type_der==COV) || (type_der==CON)) ;
	int val_res = so.val_t->get_valence() + 1 ;

	bool donames = (ind_der==' ') ? false : true ;
	bool need_sum = false ;

	Array<int> type_ind (val_res) ;
	type_ind.set(0) = COV ;
	for (int i=1 ; i<val_res ; i++)
	  type_ind.set(i) = so.val_t->get_index_type(i-1) ;

	if (donames) {
	  if (so.val_t->get_valence()>0)
		  assert (so.val_t->is_name_affected()) ;
	  // Need for summation ?
	  for (int i=1 ; i<val_res ; i++)
		  if (ind_der== so.val_t->get_name_ind()[i-1])
			  need_sum = true ;
	}


	// Flat grad space
	Term_eq fgrad (flat_grad_spher(so)) ;
	bool doder = (fgrad.der_t==0x0) ? false : true ;

	// Connexions parts 
	Base_tensor basis (so.val_t->get_space(), SPHERICAL_BASIS) ;
	Tensor auxi_val (so.val_t->get_space(), val_res, type_ind, basis) ;
	auxi_val = 0 ;
	auxi_val.std_base() ;
	
	if (donames) {
	  // Set the names of the indices :
	  auxi_val.set_name_affected() ;
	  auxi_val.set_name_ind(0, ind_der) ;
	  fgrad.val_t->set_name_affected() ;
	  fgrad.val_t->set_name_ind(0, ind_der) ;
	  for (int i=1 ; i<val_res ; i++) {
		  auxi_val.set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;
		  fgrad.val_t->set_name_ind (i, so.val_t->get_name_ind()[i-1]) ;
	  }
	}

	// Part derivative
	//Loop on the components :
	Index pos_auxi_bis(auxi_val) ;
	Index pos_so_bis (*so.val_t) ;

	// Loop indice summation on connection symbols 
	for (int ind_sum=0 ; ind_sum<val_res-1 ; ind_sum++) {

	  //Loop on the components :
	  Index pos_auxi(auxi_val) ;
	  Index pos_so (*so.val_t) ;

	  do {
		  for (int i=0 ; i<val_res-1 ; i++)
			  pos_so.set(i) = pos_auxi(i+1) ;
		  // Different cases of the derivative index :
		  switch (pos_auxi(0)) {
		      case 0 : 
			// Dr nothing
			break ;
		      case 1 :
			// Dtheta
			// Different cases of the source index
			switch (pos_auxi(ind_sum+1)) {
			    case 0 :
			      //Dtheta S_r 
			      pos_so.set(ind_sum) = 1 ;
			      auxi_val.set(pos_auxi).set_domain(num_dom) -= (*so.val_t)(pos_so)(num_dom) ;
			      break ;
			    case 1 :
			      // Dtheta S_theta
			      pos_so.set(ind_sum) = 0 ;
			      auxi_val.set(pos_auxi).set_domain(num_dom) += (*so.val_t)(pos_so)(num_dom) ;
			      break ;
			    case 2 :
			      //Dtheta S_phi 
			      break ;
			    default :
			      cerr << "Bad indice in Domain_shell_outer_adapted::derive_flat_spher" << endl  ;
			      abort() ;
			}
			break ;
		      case 2 :
			// Dphi
			// Different cases of the source index
			switch (pos_auxi(ind_sum+1)) {
			    case 0 :
			      //Dphi S_r 
			      pos_so.set(ind_sum) = 2 ;
			      auxi_val.set(pos_auxi).set_domain(num_dom) -= (*so.val_t)(pos_so)(num_dom) ;
			      break ;
			    case 1 :
			      // Dphi S_theta
			      pos_so.set(ind_sum) = 2 ;
			      auxi_val.set(pos_auxi).set_domain(num_dom) -= (*so.val_t)(pos_so)(num_dom).mult_cos_theta().div_sin_theta() ;
			      break ;
			    case 2 :
			      //Dphi S_phi
			      pos_so.set(ind_sum) = 0 ;
			      auxi_val.set(pos_auxi).set_domain(num_dom) += (*so.val_t)(pos_so)(num_dom) ;
			      pos_so.set(ind_sum) = 1 ;
			      auxi_val.set(pos_auxi).set_domain(num_dom) += (*so.val_t)(pos_so)(num_dom).mult_cos_theta().div_sin_theta() ;
			      break ;
			    default :
			      cerr << "Bad indice in Domain_shell_outer_adapted::derive_flat_spher" << endl  ;
			      abort() ;
			}
			break ;
		      default :
			cerr << "Bad indice in Domain_shell_outer_adapted::derive_flat_spher" << endl  ;
			abort() ;
		  }
	  }
	  while (pos_auxi.inc()) ;
	}
  
	if (!doder) {
		// No need for derivative :
		Term_eq occi (num_dom, auxi_val) ;
		
		Term_eq auxi (occi / (*rad_term_eq) + fgrad) ;
		
		// If derive contravariant : manipulate first indice :
		if (type_der==CON) 
			manipulator->manipulate_ind (auxi, 0) ;
	
		if (!need_sum)
			return auxi ;
		else
			return Term_eq (num_dom, auxi.val_t->do_summation_one_dom(num_dom)) ;
	}
	else {
		// Need to compute the derivative :
		// Tensor for der
		Tensor auxi_der (so.val_t->get_space(), val_res, type_ind, basis) ;
		auxi_der = 0 ;	

		if (donames) {
		  // Set the names of the indices :
		  auxi_der.set_name_affected() ;
		  auxi_der.set_name_ind(0, ind_der) ;
		  fgrad.der_t->set_name_affected() ;
		  fgrad.der_t->set_name_ind(0, ind_der) ;
		  for (int i=1 ; i<val_res ; i++) {
			  auxi_der.set_name_ind(i, so.der_t->get_name_ind()[i-1]) ;
			  fgrad.der_t->set_name_ind(i, so.der_t->get_name_ind()[i-1]) ;
		  }
		}

		// Part derivative
		//Loop on the components :
		Index pos_auxi_der_bis(auxi_der) ;

		// Loop indice summation on connection symbols 
	for (int ind_sum=0 ; ind_sum<val_res-1 ; ind_sum++) {

	  //Loop on the components :
	  Index pos_auxi_der(auxi_val) ;
	  Index pos_so (*so.val_t) ;

	  do {
		  for (int i=0 ; i<val_res-1 ; i++)
			  pos_so.set(i) = pos_auxi_der(i+1) ;
		  // Different cases of the derivative index :
		  switch (pos_auxi_der(0)) {
		      case 0 : 
			// Dr nothing
			break ;
		      case 1 :
			// Dtheta
			// Different cases of the source index
			switch (pos_auxi_der(ind_sum+1)) {
			    case 0 :
			      //Dtheta S_r 
			      pos_so.set(ind_sum) = 1 ;
			      auxi_der.set(pos_auxi_der).set_domain(num_dom) -= (*so.der_t)(pos_so)(num_dom) ;
			      break ;
			    case 1 :
			      // Dtheta S_theta
			      pos_so.set(ind_sum) = 0 ;
			      auxi_der.set(pos_auxi_der).set_domain(num_dom) += (*so.der_t)(pos_so)(num_dom) ;
			      break ;
			    case 2 :
			      //Dtheta S_phi 
			      break ;
			    default :
			      cerr << "Bad indice in Domain_shell_outer_adapted::derive_flat_spher" << endl  ;
			      abort() ;
			}
			break ;
		      case 2 :
			// Dphi
			// Different cases of the source index
			switch (pos_auxi_der(ind_sum+1)) {
			    case 0 :
			      //Dphi S_r 
			      pos_so.set(ind_sum) = 2 ;
			      auxi_der.set(pos_auxi_der).set_domain(num_dom) -= (*so.der_t)(pos_so)(num_dom) ;
			      break ;
			    case 1 :
			      // Dphi S_theta
			      pos_so.set(ind_sum) = 2 ;
			      auxi_der.set(pos_auxi_der).set_domain(num_dom) -= (*so.der_t)(pos_so)(num_dom).mult_cos_theta().div_sin_theta() ;
			      break ;
			    case 2 :
			      //Dphi S_phi
			      pos_so.set(ind_sum) = 0 ;
			      auxi_der.set(pos_auxi_der).set_domain(num_dom) += (*so.der_t)(pos_so)(num_dom) ;
			      pos_so.set(ind_sum) = 1 ;
			      auxi_der.set(pos_auxi_der).set_domain(num_dom) += (*so.der_t)(pos_so)(num_dom).mult_cos_theta().div_sin_theta() ;
			      break ;
			    default :
			      cerr << "Bad indice in Domain_shell_outer_adapted::derive_flat_spher" << endl  ;
			      abort() ;
			}
			break ;
		      default :
			cerr << "Bad indice in Domain_shell_outer_adapted::derive_flat_spher" << endl  ;
			abort() ;
		  }
	  }
	  while (pos_auxi_der.inc()) ;
	}

		// Need for derivative :
		Term_eq occi (num_dom, auxi_val, auxi_der) ;	
		Term_eq auxi (occi/(*rad_term_eq) + fgrad) ;

		


		// If derive contravariant : manipulate first indice :
		if (type_der==CON)
			manipulator->manipulate_ind (auxi, 0) ;
		if (!need_sum)
			return auxi ;
		else
			return Term_eq (num_dom, auxi.val_t->do_summation_one_dom(num_dom), 
							auxi.der_t->do_summation_one_dom(num_dom)) ;
	}
}


Term_eq Domain_shell_outer_adapted::derive_flat_cart (int type_der, char ind_der, const Term_eq& so, const Metric* manipulator) const {
 
	bool doder = (so.der_t==0x0) ? false : true ;
	if ((type_der==CON) && (manipulator->p_met_con[num_dom]->der_t==0x0))
	  doder = false ;

	assert ((type_der==COV) || (type_der==CON)) ;
	int val_res = so.val_t->get_valence() + 1 ;

	bool doname = true ;
	if (so.val_t->get_valence()>0)
		if (!so.val_t->is_name_affected()) 
			doname = false;

	Term_eq fgrad (partial_cart(so)) ;
	if (fgrad.der_t==0x0)
	  doder = false ;
	
	// Need for summation ?
	bool need_sum = false ;
	if (doname)
		for (int i=1 ; i<val_res ; i++)
			if (ind_der== so.val_t->get_name_ind()[i-1])
				need_sum = true ;
	
	// Set the names of the indices :
	if (doname) {
		fgrad.val_t->set_name_affected() ;
		fgrad.val_t->set_name_ind(0, ind_der) ;
		for (int i=1 ; i<val_res ; i++)
			fgrad.val_t->set_name_ind(i, so.val_t->get_name_ind()[i-1]) ;
	}

	if (!doder) {
		// If derive contravariant : manipulate first indice :
		if (type_der==CON) 
			manipulator->manipulate_ind (fgrad, 0) ;
	
		if (!need_sum)
			return fgrad ;
		else
			return Term_eq (num_dom, fgrad.val_t->do_summation_one_dom(num_dom)) ;
	}
	else {
		// Need to compute the derivative :
		
		fgrad.der_t->set_name_affected() ;
		fgrad.der_t->set_name_ind(0, ind_der) ;
		for (int i=1 ; i<val_res ; i++)
			fgrad.der_t->set_name_ind(i, so.der_t->get_name_ind()[i-1]) ;

		// If derive contravariant : manipulate first indice :
		if (type_der==CON)
		    manipulator->manipulate_ind (fgrad, 0) ;
		
		if (!need_sum)
			return fgrad ;
		else 
		  return Term_eq (num_dom, fgrad.val_t->do_summation_one_dom(num_dom), 
							fgrad.der_t->do_summation_one_dom(num_dom)) ;
	}
 
}
}

