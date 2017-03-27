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
#include "adapted_polar.hpp"
#include "point.hpp"
#include "array_math.cpp"
#include "val_domain.hpp"
#include "scalar.hpp"
#include "vector.hpp"

namespace Kadath {
Term_eq Domain_polar_shell_outer_adapted::dr_term_eq (const Term_eq& so) const {
  return derive_r(so) ;
}
  
Term_eq Domain_polar_shell_outer_adapted::derive_r (const Term_eq& so) const {
  
  assert (so.get_dom() == num_dom) ;
  
  Tensor val_res(*so.val_t, false) ;
  
  for (int cmp = 0 ; cmp<so.val_t->get_n_comp() ; cmp++) {
    Array<int>ind (val_res.indices(cmp)) ;
    val_res.set(ind).set_domain(num_dom) = (*so.val_t)(ind)(num_dom).der_var(1) / (*der_rad_term_eq->val_t)()(num_dom)  ;
  }
  
  bool doder = ((so.der_t==0x0) || (der_rad_term_eq->der_t==0x0)) ? false : true ;
  if (doder) {
      Tensor der_res (*so.der_t, false) ;
    for (int cmp = 0 ; cmp<so.val_t->get_n_comp() ; cmp++) {
      Array<int>ind (val_res.indices(cmp)) ;
      der_res.set(ind).set_domain(num_dom) = ((*so.der_t)(ind)(num_dom).der_var(1)*(*der_rad_term_eq->val_t)()(num_dom) 
		  - (*so.val_t)(ind)(num_dom).der_var(1)*(*der_rad_term_eq->der_t)()(num_dom)) / 
		  ((*der_rad_term_eq->val_t)()(num_dom) * (*der_rad_term_eq->val_t)()(num_dom))  ; 
    }
    return Term_eq (num_dom, val_res, der_res) ;
  }
  else
    return Term_eq (num_dom, val_res) ;
}

Term_eq Domain_polar_shell_outer_adapted::derive_t (const Term_eq& so) const {

   assert (so.get_dom() == num_dom) ;
  Term_eq* dtprime ;
  Tensor val_res(*so.val_t, false) ;
  
  for (int cmp = 0 ; cmp<so.val_t->get_n_comp() ; cmp++) {
    Array<int>ind (val_res.indices(cmp)) ;
    val_res.set(ind).set_domain(num_dom) = (*so.val_t)(ind)(num_dom).der_var(2)  ;
  }
  bool doder = (so.der_t==0x0) ? false : true ;
 
  if (doder) {
    Tensor der_res (*so.der_t, false) ;
     for (int cmp = 0 ; cmp<so.val_t->get_n_comp() ; cmp++) {
       Array<int>ind (val_res.indices(cmp)) ;
       der_res.set(ind).set_domain(num_dom) = (*so.der_t)(ind)(num_dom).der_var(2)  ; 
     }
    dtprime = new Term_eq (num_dom, val_res, der_res) ;
  }
  else
    dtprime = new Term_eq (num_dom, val_res) ;
  
  
  Term_eq res (*dtprime - (*dt_rad_term_eq) * derive_r(so)) ;
  delete dtprime ;
  return res ;
  
}

Term_eq Domain_polar_shell_outer_adapted::flat_grad_spher (const Term_eq& so) const {
   
  assert (so.get_dom() == num_dom) ;
  
  int valso = so.get_val_t().get_valence() ;
  Array<int> type_ind (valso+1) ; 
  type_ind.set(0) = COV ;
  for (int i=0 ; i<valso ; i++)
    type_ind.set(i+1) = so.get_val_t().get_index_type(i) ;
  
  Base_tensor basis (sp) ;
  basis.set_basis(num_dom) =  SPHERICAL_BASIS ;
  Tensor res (sp, valso+1 , type_ind, basis) ;
  
  Term_eq comp_r (derive_r(so)) ;
  Term_eq comp_t (derive_t(so)/(*rad_term_eq)) ;
  
  // Loop on cmp :
  for (int nc=0 ; nc<so.get_val_t().get_n_comp() ; nc++) {
    
    Array<int> ind (so.get_val_t().indices(nc)) ;
    Array<int> indtarget (valso+1) ;
    for (int i=0 ; i<valso ; i++)
	indtarget.set(i+1) = ind(i) ;
    
    // R comp :
    indtarget.set(0) = 1 ;
    res.set(indtarget).set_domain(num_dom) = comp_r.get_val_t()(ind)(num_dom);
    // theta comp :
    indtarget.set(0) = 2 ;
    res.set(indtarget).set_domain(num_dom) = comp_t.get_val_t()(ind)(num_dom) ;
  }
  bool doder = ((so.der_t==0x0) || (outer_radius_term_eq->der_t==0x0)) ? false : true ;
  
  if (!doder) {
    return Term_eq (num_dom, res) ;
  }
  else {
    
    Tensor resder (sp, valso+1 , type_ind, basis) ;
     // Loop on cmp :
  for (int nc=0 ; nc<so.get_val_t().get_n_comp() ; nc++) {
    
    Array<int> ind (so.get_val_t().indices(nc)) ;
    Array<int> indtarget (valso+1) ;
    for (int i=0 ; i<valso ; i++)
	indtarget.set(i+1) = ind(i) ;
   
    // R comp :
    indtarget.set(0) = 1 ;
    resder.set(indtarget).set_domain(num_dom) = comp_r.get_der_t()(ind)(num_dom);
    // theta comp :
    indtarget.set(0) = 2 ;
    resder.set(indtarget).set_domain(num_dom) = comp_t.get_der_t()(ind)(num_dom) ;
    }
   return Term_eq (num_dom, res, resder) ; 
  }
}

void Domain_polar_shell_outer_adapted::do_normal_spher () const {
  
    Term_eq grad (flat_grad_spher(*rad_term_eq - *outer_radius_term_eq)) ;

    Scalar val_norme (sp) ;
    val_norme.set_domain(num_dom) = sqrt((*grad.val_t)(1)(num_dom)*(*grad.val_t)(1)(num_dom)
			      + (*grad.val_t)(2)(num_dom)*(*grad.val_t)(2)(num_dom)) ;
    bool doder = ((grad.der_t==0x0) || (outer_radius_term_eq->der_t==0x0)) ? false : true ;
    if (doder) {
       Scalar der_norme (sp) ;
       der_norme.set_domain(num_dom) = ((*grad.der_t)(1)(num_dom)*(*grad.val_t)(1)(num_dom)
			      + (*grad.der_t)(2)(num_dom)*(*grad.val_t)(2)(num_dom)) / val_norme(num_dom) ;			      
	Term_eq norme (num_dom, val_norme, der_norme) ;
	if (normal_spher==0x0)
	  normal_spher = new Term_eq (grad / norme) ;
	else
	  *normal_spher = Term_eq(grad/norme) ;
    }
    else {
      Term_eq norme (num_dom, val_norme) ;
      if (normal_spher==0x0)
	normal_spher = new Term_eq (grad / norme)  ;
      else
	*normal_spher = Term_eq(grad/norme) ;
    }
}

void Domain_polar_shell_outer_adapted::do_normal_cart () const {
  
  
    do_normal_spher() ;
    
    Base_tensor basis (sp) ;
    basis.set_basis(num_dom) =  CARTESIAN_BASIS ;
    Vector val (sp, CON, basis) ;
    
    val.set(1).set_domain(num_dom) = mult_sin_theta((*normal_spher->val_t)(1)(num_dom)) + mult_cos_theta((*normal_spher->val_t)(2)(num_dom)) ;
    val.set(2).set_domain(num_dom) = mult_cos_theta((*normal_spher->val_t)(1)(num_dom)) - mult_sin_theta((*normal_spher->val_t)(2)(num_dom)) ;
    
    bool doder = (normal_spher->der_t==0x0) ? false : true ;
    if (doder) {
       Vector der (sp, CON, basis) ;
    
    der.set(1).set_domain(num_dom) = mult_sin_theta((*normal_spher->der_t)(1)(num_dom)) + mult_cos_theta((*normal_spher->der_t)(2)(num_dom)) ;
    der.set(2).set_domain(num_dom) = mult_cos_theta ((*normal_spher->der_t)(1)(num_dom)) - mult_sin_theta((*normal_spher->der_t)(2)(num_dom)) ;
      
	if (normal_cart==0x0)
	  normal_cart = new Term_eq (num_dom, val, der) ;
	else
	  *normal_cart = Term_eq(num_dom, val,der) ;
    }
    else {
     if (normal_cart==0x0)
	  normal_cart = new Term_eq (num_dom, val) ;
	else
	  *normal_cart = Term_eq(num_dom, val) ;
    }
}

Term_eq Domain_polar_shell_outer_adapted::der_normal_term_eq (const Term_eq& so, int bound) const {
  
  
  switch (bound) {
    case OUTER_BC : {    
     // Deformed surface
      if (normal_spher == 0x0)
	do_normal_spher() ;
      
      int valso = so.get_val_t().get_valence() ;
      
      Term_eq grad (flat_grad_spher(so)) ;
  
      Tensor res (so.get_val_t(), false) ;
      
      // Loop on cmp :
      for (int nc=0 ; nc<so.get_val_t().get_n_comp() ; nc++) {
    
	Array<int> ind (so.get_val_t().indices(nc)) ;
	Array<int> indgrad (valso+1) ;
	for (int i=0 ; i<valso ; i++)
	    indgrad.set(i+1) = ind(i) ;
	
	indgrad.set(0) = 1 ;
	res.set(ind).set_domain(num_dom) = (*grad.val_t)(indgrad)(num_dom) * (*normal_spher->val_t)(1)(num_dom) ;
	indgrad.set(0) = 2 ;
	res.set(ind).set_domain(num_dom) += (*grad.val_t)(indgrad)(num_dom) * (*normal_spher->val_t)(2)(num_dom) ;
      }
			  
      bool doder = ((grad.der_t == 0x0) || (normal_spher->der_t == 0x0)) ? false : true ;
      if (doder) {
	
	  Tensor der (so.get_val_t(), false) ;
      
      // Loop on cmp :
      for (int nc=0 ; nc<so.get_val_t().get_n_comp() ; nc++) {
    
	Array<int> ind (so.get_val_t().indices(nc)) ;
	Array<int> indgrad (valso+1) ;
	for (int i=0 ; i<valso ; i++)
	    indgrad.set(i+1) = ind(i) ;
	
	indgrad.set(0) = 1 ;
	der.set(ind).set_domain(num_dom) = (*grad.der_t)(indgrad)(num_dom) * (*normal_spher->val_t)(1)(num_dom) 
					+ (*grad.val_t)(indgrad)(num_dom) * (*normal_spher->der_t)(1)(num_dom) ;
	indgrad.set(0) = 2 ;
	der.set(ind).set_domain(num_dom) += (*grad.der_t)(indgrad)(num_dom) * (*normal_spher->val_t)(2)(num_dom) 
					+ (*grad.val_t)(indgrad)(num_dom) * (*normal_spher->der_t)(2)(num_dom) ;
      }
			  
	return Term_eq (num_dom, res, der) ;
      }
      else 
	return Term_eq (num_dom, res) ;
    }
    case INNER_BC : {
	return derive_r(so) ;
    }
    default : 
	cerr << "Unknown boundary in Domain_polar_shell_outer_adapted::der_normal" << endl ;
	abort() ;
  }
}


Term_eq Domain_polar_shell_outer_adapted::lap_term_eq (const Term_eq& so, int mm) const {
   
   if (so.get_val_t().get_valence() != 0) {
      cerr << "Domain_polar_shell_outer_adapted::lap_term_eq only defined for scalars" << endl ;
      abort() ;
  }
  
  // Angular part :
  
  Term_eq auxi (do_comp_by_comp(so, &Domain::div_sin_theta)) ;
 
  Term_eq dert (derive_t(so));
  Term_eq p1 (derive_t(dert)) ;
  Term_eq p2 (do_comp_by_comp(do_comp_by_comp(dert, &Domain::mult_cos_theta) -mm* auxi, &Domain::div_sin_theta)) ;

  Term_eq dr(derive_r(so)) ;
  
  Term_eq res (derive_r(dr) + 2 *dr / (*rad_term_eq) + (p1+p2)/(*rad_term_eq)/(*rad_term_eq));    
 
  
  return res ;
}

Term_eq Domain_polar_shell_outer_adapted::lap2_term_eq (const Term_eq& so, int mm) const {
    if (mm!=0) {
      cerr << "Domain_polar_shell_outer_adapted::lap2_term_eq not defined for m != 0 (for now)" << endl ;
      abort() ;
  }
   
   if (so.get_val_t().get_valence() != 0) {
      cerr << "Domain_polar_shell_outer_adapted::lap2_term_eq only defined for scalars" << endl ;
      abort() ;
  }
  
  // Angular part :
  Term_eq dert (derive_t(so));
  Term_eq p1 (derive_t(dert)) ;
  
  Term_eq dr(derive_r(so)) ;
  
  Term_eq res (derive_r(dr) + dr / (*rad_term_eq) + p1/(*rad_term_eq)/(*rad_term_eq));    
 
  
  return res ;
}

Term_eq Domain_polar_shell_outer_adapted::mult_r_term_eq (const Term_eq& so) const {
 
  
  return so * (*rad_term_eq) ;
}

Term_eq Domain_polar_shell_outer_adapted::div_r_term_eq (const Term_eq& so) const {
 
  
  return so / (*rad_term_eq) ;
}

void Domain_polar_shell_outer_adapted::update_term_eq (Term_eq* so) const {
  
  Tensor der (*so->val_t, false) ;
  for (int cmp=0 ; cmp<so->val_t->get_n_comp() ; cmp ++) {
    
    
    
    Val_domain derr ((*(*so->val_t).cmp[cmp])(num_dom).der_var(1) / (*der_rad_term_eq->val_t)()(num_dom)) ;
    
    if (!derr.check_if_zero()) { 
      Val_domain res(this) ;
      res.allocate_conf() ;
      Index pos (nbr_points) ;     
      do {
	res.set(pos) =  (*(*outer_radius_term_eq->der_t).cmp[0])(num_dom)(pos) / 2. * (1+((*coloc[0])(pos(0)))) * derr(pos) ;
      }
    while (pos.inc()) ;
  res.set_base() = (*(*so->val_t).cmp[cmp])(num_dom).get_base() ;
  der.cmp[cmp]->set_domain(num_dom) = res ;
    }
    else
      der.cmp[cmp]->set_domain(num_dom).set_zero() ;
  }
  
  so->set_der_t (der) ;
}

Term_eq Domain_polar_shell_outer_adapted::partial_spher (const Term_eq& so) const {
   int dom = so.get_dom() ;
  assert (dom == num_dom) ;
  
  int valence = so.val_t->get_valence() ;

  Array<int> type_ind (valence+1) ;
  type_ind.set(0) = COV ;
  for (int i=1 ; i<valence+1 ; i++)
      type_ind.set(i) = so.val_t->get_index_type(i-1) ;
  
  Base_tensor basis (so.get_val_t().get_space()) ;
  basis.set_basis(dom) = SPHERICAL_BASIS ;

  Term_eq comp_r (derive_r(so)) ;
  Term_eq comp_t (derive_t(so) / (*rad_term_eq)) ;
  
  Tensor val_res (so.get_val_t().get_space(), valence+1, type_ind, basis) ;
  {
  Index pos_so (*so.val_t) ;
  Index pos_res (val_res) ;
  do {
    for (int i=1 ; i<valence+1 ; i++)
      pos_res.set(i) = pos_so(i-1) ;
    // R part
    pos_res.set(0) = 0 ;
    val_res.set(pos_res).set_domain(num_dom) = (*comp_r.val_t)(pos_so)(num_dom) ;
    // Theta part
    pos_res.set(0) = 1 ;
    val_res.set(pos_res).set_domain(num_dom) = (*comp_t.val_t)(pos_so)(num_dom) ;
  }
  while (pos_so.inc()) ;
}
  
  bool doder = ((so.der_t==0x0) || (outer_radius_term_eq->der_t==0x0)) ? false : true ;
  if (doder) {
     Tensor der_res (so.get_val_t().get_space(), valence+1, type_ind, basis) ;
     Index pos_so (*so.der_t) ;
     Index pos_res (val_res) ;
  do {
    for (int i=1 ; i<valence+1 ; i++)
      pos_res.set(i) = pos_so(i-1) ;
    // R part
    pos_res.set(0) = 0 ;
    der_res.set(pos_res).set_domain(num_dom) = (*comp_r.der_t)(pos_so)(num_dom) ;
    // Theta part
    pos_res.set(0) = 1 ;
    der_res.set(pos_res).set_domain(num_dom) = (*comp_t.der_t)(pos_so)(num_dom) ;
  }
  while (pos_so.inc()) ;
    
    return Term_eq (num_dom, val_res, der_res) ;
  }
  else
      return Term_eq (num_dom, val_res) ;
}


Term_eq Domain_polar_shell_outer_adapted::partial_cart (const Term_eq& so) const {
  int dom = so.get_dom() ;
  assert (dom == num_dom) ;
  
  int valence = so.val_t->get_valence() ;

  Array<int> type_ind (valence+1) ;
  type_ind.set(0) = COV ;
  for (int i=1 ; i<valence+1 ; i++)
      type_ind.set(i) = so.val_t->get_index_type(i-1) ;
  
  Base_tensor basis (so.get_val_t().get_space()) ;
  basis.set_basis(dom) = CARTESIAN_BASIS ;

  Term_eq comp_r (derive_r(so)) ;
  Term_eq comp_t (derive_t(so) / (*rad_term_eq)) ;
  
  Tensor val_res (so.get_val_t().get_space(), valence+1, type_ind, basis) ;
  {
  Index pos_so (*so.val_t) ;
  Index pos_res (val_res) ;
  do {
    for (int i=1 ; i<valence+1 ; i++)
      pos_res.set(i) = pos_so(i-1) ;
   
    // X part
    pos_res.set(0) = 0 ;
    val_res.set(pos_res).set_domain(num_dom) = mult_sin_theta((*comp_r.val_t)(pos_so)(num_dom)) + mult_cos_theta((*comp_t.val_t)(pos_so)(num_dom)) ;
    // Z part
    pos_res.set(0) = 1 ;
    val_res.set(pos_res).set_domain(num_dom) = 
      mult_cos_theta((*comp_r.val_t)(pos_so)(num_dom)) - mult_sin_theta((*comp_t.val_t)(pos_so)(num_dom)) ;
  }
  while (pos_so.inc()) ;
}
  
  bool doder = ((so.der_t==0x0) || (rad_term_eq->der_t==0x0)) ? false : true ;
  if (doder) {
     Tensor der_res (so.get_val_t().get_space(), valence+1, type_ind, basis) ;
     Index pos_so (*so.der_t) ;
     Index pos_res (val_res) ;
  do {
    for (int i=1 ; i<valence+1 ; i++)
      pos_res.set(i) = pos_so(i-1) ;
    
    Val_domain auxi (mult_sin_theta((*comp_r.der_t)(pos_so)(num_dom)) + mult_cos_theta((*comp_t.der_t)(pos_so)(num_dom))) ;
    // X part
    pos_res.set(0) = 0 ;
    der_res.set(pos_res).set_domain(num_dom) =  mult_sin_theta((*comp_r.der_t)(pos_so)(num_dom)) + mult_cos_theta((*comp_t.der_t)(pos_so)(num_dom)) ;
    // Z part
    pos_res.set(0) = 1 ;
    der_res.set(pos_res).set_domain(num_dom) = 
      mult_cos_theta((*comp_r.der_t)(pos_so)(num_dom)) - mult_sin_theta((*comp_t.der_t)(pos_so)(num_dom)) ;
  }
  while (pos_so.inc()) ;
    
    return Term_eq (num_dom, val_res, der_res) ;
  }
  else
      return Term_eq (num_dom, val_res) ;
}

Term_eq Domain_polar_shell_outer_adapted::grad_term_eq (const Term_eq& so) const {
  return partial_cart(so) ;
}

const Term_eq* Domain_polar_shell_outer_adapted::give_normal (int bound, int tipe) const {
  assert (bound==INNER_BC) ;
  switch (tipe) {
    case CARTESIAN_BASIS : 
      if (normal_cart==0x0)
	do_normal_cart() ;
      return normal_cart ;
    
    case SPHERICAL_BASIS : 
      if (normal_spher==0x0)
	do_normal_spher() ;
      return normal_spher ;
      
    default: 
      cerr << "Unknown type of tensorial basis in Domain_polar_shell_outer_adapted::give_normal" << endl ;
      abort() ;
  }
}

Term_eq Domain_polar_shell_outer_adapted::integ_term_eq (const Term_eq&, int) const {
  cerr << "Domain_polar_shell_outer_adapted::integ_term_eq not implemented yet" << endl ;
  abort() ;
}
}

