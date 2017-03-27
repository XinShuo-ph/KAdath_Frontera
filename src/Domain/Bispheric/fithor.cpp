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

#include "bispheric.hpp"
#include "param.hpp"
#include "term_eq.hpp"
#include "scalar.hpp"
#include "vector.hpp"

namespace Kadath {
Val_domain Domain_bispheric_rect::fithor (const Val_domain& so) const {
  
  Val_domain constante (this) ;
  constante.allocate_conf() ;
  
  Index pos (nbr_points) ;
  Index poshor (nbr_points) ;
  do {
    poshor.set(0) = 0 ;
    poshor.set(1) = pos(1) ;
    poshor.set(2) = pos(2) ;
    
    constante.set(pos) = so(poshor) ;
  }
  while (pos.inc()) ;
  
  double rhor = aa / sinh(fabs(eta_minus)) ;
  double xc = aa*cosh(eta_minus)/sinh(eta_minus) ;
  
  Val_domain radloc2 ((get_cart(1)-xc)*(get_cart(1)-xc)+get_cart(2)*get_cart(2)+get_cart(3)*get_cart(3)) ;
   
  Val_domain result (constante*rhor*rhor/radloc2) ;
  result.set_base() = so.get_base() ; 
  return result ;
}


Term_eq Domain_bispheric_rect::fithor (const Term_eq& so) const {
    // Check it is a tensor
    if (so.get_type_data() != TERM_T) {
		cerr << "fithor only defined with respect for a tensor" << endl ;
		abort() ;
    }
  
    // Right domain ?
    int dom = so.get_dom() ;
    if (this != so.get_val_t().get_space().get_domain(dom)) {
	cerr << "Domain mismatch in Domain_bispheric_rect::fithor (Term_eq version)" << endl ;
	abort() ;
    }
   int valence = so.get_val_t().get_valence() ;
  
    // Scalar ?
    if ((valence!=0) && (valence!=1))  {
	cerr << "Domain_bispheric_chi_first::fithor not definded for valence= " << valence << endl ;
	abort() ;
    }

    if (valence==0) {
      Scalar resval (so.get_val_t(), false) ;
      Val_domain value (so.get_val_t()()(dom)) ;
      if (value.check_if_zero()) 
	resval.set().set_domain(dom).set_zero() ;
      else
	  resval.set().set_domain(dom) = fithor(value) ;	
    
      if (so.get_p_der_t() !=0x0) {
	Scalar resder (so.get_der_t(), false) ;
	Val_domain valder (so.get_der_t()()(dom)) ;
	if (valder.check_if_zero())
	  resder.set().set_domain(dom).set_zero() ;
	else 
	  resder.set().set_domain(dom) = fithor(valder) ;
	return Term_eq (dom, resval, resder) ;
      }
      else return Term_eq (dom, resval) ;
    }
    
    
    if (valence==1) {
        Vector resval (so.get_val_t().get_space(), so.get_val_t().get_index_type(0), so.get_val_t().get_basis()) ;
	for (int i=1 ; i<=3 ; i++) {
	    Val_domain value (so.get_val_t()(i)(dom)) ;
	    if (value.check_if_zero()) 
	      resval.set(i).set_domain(dom).set_zero() ;
	    else
	      resval.set(i).set_domain(dom) = fithor(value) ;	
	}
	
	if (so.get_val_t().is_name_affected()) {
	    resval.set_name_affected() ;
	    resval.set_name_ind(0, so.get_val_t().get_name_ind()[0]) ;
	}
	
	if (so.get_p_der_t() !=0x0) {
	  Vector resder (so.get_der_t().get_space(), so.get_der_t().get_index_type(0), so.get_der_t().get_basis()) ;
	for (int i=1 ; i<=3 ; i++) {
	  Val_domain valder (so.get_der_t()(i)(dom)) ;
	  if (valder.check_if_zero())
	    resder.set(i).set_domain(dom).set_zero() ;
	  else 
	    resder.set(i).set_domain(dom) = fithor(valder) ;
	} 
	if (so.get_der_t().is_name_affected()) {
	    resder.set_name_affected() ;
	    resder.set_name_ind(0, so.get_val_t().get_name_ind()[0]) ;
	}
	  return Term_eq (dom, resval, resder) ;
      }
      else return Term_eq (dom, resval) ;
    }
    
    abort() ;
}

Val_domain Domain_bispheric_chi_first::fithor (const Val_domain& so) const {
  
 
  Val_domain constante (this) ;
  constante.allocate_conf() ;
  
  Index pos (nbr_points) ;
  Index poshor (nbr_points) ;
  do {
    poshor.set(0) = 0 ;
    poshor.set(1) = pos(1) ;
    poshor.set(2) = pos(2) ;
    
    constante.set(pos) = so(poshor) ;
  }
  while (pos.inc()) ;
 
  double rhor = aa / sinh(fabs(eta_lim)) ;
  double xc = aa*cosh(eta_lim)/sinh(eta_lim) ;
   
  Val_domain radloc2 ((get_cart(1)-xc)*(get_cart(1)-xc)+get_cart(2)*get_cart(2)+get_cart(3)*get_cart(3)) ;
   
  
  Val_domain result (constante*rhor*rhor/radloc2) ;
  result.set_base() = so.get_base() ;
  return result ;
}


Term_eq Domain_bispheric_chi_first::fithor (const Term_eq& so) const {
    // Check it is a tensor
    if (so.get_type_data() != TERM_T) {
		cerr << "fithor only defined with respect for a tensor" << endl ;
		abort() ;
    }
  
    // Right domain ?
    int dom = so.get_dom() ;
    if (this != so.get_val_t().get_space().get_domain(dom)) {
	cerr << "Domain mismatch in Domain_bispheric_chi_first::fithor (Term_eq version)" << endl ;
	abort() ;
    }

    int valence = so.get_val_t().get_valence() ;
  
    // Scalar ?
    if ((valence!=0) && (valence!=1))  {
	cerr << "Domain_bispheric_chi_first::fithor not definded for valence= " << valence << endl ;
	abort() ;
    }

    if (valence==0) {
      Scalar resval (so.get_val_t(), false) ;
      Val_domain value (so.get_val_t()()(dom)) ;
      if (value.check_if_zero()) 
	resval.set().set_domain(dom).set_zero() ;
      else
	  resval.set().set_domain(dom) = fithor(value) ;	
    
      if (so.get_p_der_t() !=0x0) {
	Scalar resder (so.get_der_t(), false) ;
	Val_domain valder (so.get_der_t()()(dom)) ;
	if (valder.check_if_zero())
	  resder.set().set_domain(dom).set_zero() ;
	else 
	  resder.set().set_domain(dom) = fithor(valder) ;
	return Term_eq (dom, resval, resder) ;
      }
      else return Term_eq (dom, resval) ;
    }
    
    
    if (valence==1) {
        Vector resval (so.get_val_t().get_space(), so.get_val_t().get_index_type(0), so.get_val_t().get_basis()) ;
	for (int i=1 ; i<=3 ; i++) {
	    Val_domain value (so.get_val_t()(i)(dom)) ;
	    if (value.check_if_zero()) 
	      resval.set(i).set_domain(dom).set_zero() ;
	    else
	      resval.set(i).set_domain(dom) = fithor(value) ;	
	}
	if (so.get_val_t().is_name_affected()) {
	    resval.set_name_affected() ;
	    resval.set_name_ind(0, so.get_val_t().get_name_ind()[0]) ;
	}
	
	if (so.get_p_der_t() !=0x0) {
	  Vector resder (so.get_der_t().get_space(), so.get_der_t().get_index_type(0), so.get_der_t().get_basis()) ;
	for (int i=1 ; i<=3 ; i++) {
	  Val_domain valder (so.get_der_t()(i)(dom)) ;
	  if (valder.check_if_zero())
	    resder.set(i).set_domain(dom).set_zero() ;
	  else 
	    resder.set(i).set_domain(dom) = fithor(valder) ;
	}
	
	if (so.get_der_t().is_name_affected()) {
	    resder.set_name_affected() ;
	    resder.set_name_ind(0, so.get_der_t().get_name_ind()[0]) ;
	}
	  return Term_eq (dom, resval, resder) ;
      }
      else return Term_eq (dom, resval) ;
    }
    
    abort() ;
    
}
}

