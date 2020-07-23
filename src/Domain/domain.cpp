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

#include "Kadath_point_h/kadath.hpp"

//num_dom
//        ndim
//nbr_points
//        nbr_coefs
//type_base
//        coloc
//absol
//        cart
//radius
//        cart_surr
namespace Kadath {


// Constructor by copy
Domain::Domain (const Domain& so) :
        num_dom{so.num_dom}, ndim{so.ndim}, nbr_points{so.nbr_points}, nbr_coefs{so.nbr_coefs}, type_base{so.type_base},
        coloc{ndim,do_not_initialize}, absol{ndim,do_not_initialize}, cart{ndim,do_not_initialize},
        radius{so.radius ? new Val_domain{*so.radius} : nullptr}, cart_surr{ndim,do_not_initialize}
{
	for (int i=0 ; i<ndim ; i++) coloc[i] = (so.coloc[i] ? new Array<double>{*so.coloc[i]} : nullptr) ;
	for (int i=0 ; i<ndim ; i++) absol[i] = (so.absol[i] ? new Val_domain{*so.absol[i]} : nullptr) ;
	for (int i=0 ; i<ndim ; i++) cart[i] = (so.cart[i] ? new Val_domain{*so.cart[i]} : nullptr) ;
	for (int i=0 ; i<ndim ; i++) cart_surr[i] = (so.cart_surr[i] ? new Val_domain {*so.cart_surr[i]} : nullptr) ;
}


Domain::Domain (int num, FILE* fd) :
    num_dom(num), ndim{}, nbr_points(fd), nbr_coefs(fd), type_base{}, coloc{}, absol{}, cart{}, radius{nullptr},
    cart_surr{}
{
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;
	assert (ndim==nbr_points.get_ndim()) ;
	assert (ndim==nbr_coefs.get_ndim()) ;
	coloc.resize(ndim) ;
	for (auto & x : coloc) x = nullptr;
	absol.resize(ndim) ;
	for (auto & x : absol) x = nullptr ;
	cart.resize(ndim) ;
	for (auto & x : cart) x = nullptr ;
	cart_surr.resize(ndim) ;
	for (auto & x : cart_surr) x = nullptr ;
}



void Domain::save (FILE* fd) const {
	nbr_points.save(fd) ;
	nbr_coefs.save(fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
}


// destructor of the derived quantities
void Domain::del_deriv()  {
	for (int l=0 ; l<ndim ; l++) {
	    safe_delete(coloc[l]);
	    safe_delete(absol[l]);
        safe_delete(cart[l]);
        safe_delete(cart_surr[l]);
	}
	safe_delete(radius);
}

Val_domain Domain::laplacian (const Val_domain& so, int m) const {

  if (m!=0) {
    cerr << "In the general case, Laplacian must be called with m=0" << endl ;
    abort() ;
  }

  Val_domain res (so.der_abs(1).der_abs(1)) ;
  for (int j=2 ; j<=get_ndim() ; j++)
    res += so.der_abs(j).der_abs(j) ;
  return res ;
}

Val_domain Domain::laplacian2 (const Val_domain&, int) const {
  cerr << "laplacian2 not implemented for" << endl ;
  cerr << *this << endl ;
  abort() ;
 
}

Term_eq Domain::import (int numdom, int bound, int n_ope, Term_eq** parts) const {
 
  // get indices of domains :
  Array<int> zedoms(n_ope) ;
  for (int i=0 ; i<n_ope ; i++)
      zedoms.set(i) = parts[i]->get_dom() ;
  
  // Get val part :
  Tensor** parts_val = new Tensor* [n_ope] ;
  for (int i=0 ; i<n_ope ; i++)
    parts_val[i] = parts[i]->val_t ;
  
  Tensor res_val (import(numdom, bound, n_ope, zedoms, parts_val)) ;
  delete [] parts_val ;
  
  // Do the derivative parts ?
  bool doder = true ;
  for (int i=0 ; i<n_ope ; i++)
      if (parts[i]->der_t==nullptr)
	  doder = false ;
          
  if (doder) {
    Tensor** parts_der = new Tensor* [n_ope] ;
    for (int i=0 ; i<n_ope ; i++)
      parts_der[i] = parts[i]->der_t ;
    Tensor res_der (import(numdom, bound, n_ope, zedoms, parts_der)) ;
    delete [] parts_der ;
    return Term_eq (numdom, res_val, res_der) ;
  }
  else
    return Term_eq (numdom, res_val) ;
}

// Output operator
ostream& operator<< (ostream& o, const Domain& so) {
       const Domain_nucleus* nuc = dynamic_cast<const Domain_nucleus*>(&so) ;
       const Domain_shell* shell = dynamic_cast<const Domain_shell*>(&so) ;
       const Domain_shell_log* shell_log = dynamic_cast<const Domain_shell_log*>(&so) ; 
       const Domain_shell_surr* shell_surr = dynamic_cast<const Domain_shell_surr*>(&so) ;
       const Domain_compact* compact = dynamic_cast<const Domain_compact*>(&so) ;  
       const Domain_bispheric_rect* bi_rect = dynamic_cast<const Domain_bispheric_rect*>(&so) ;
       const Domain_bispheric_chi_first* bi_chi = dynamic_cast<const Domain_bispheric_chi_first*>(&so) ;       
       const Domain_bispheric_eta_first* bi_eta = dynamic_cast<const Domain_bispheric_eta_first*>(&so) ;
       const Domain_critic_inner* critic_inner = dynamic_cast<const Domain_critic_inner*>(&so) ;       
       const Domain_critic_outer* critic_outer = dynamic_cast<const Domain_critic_outer*>(&so) ;
       const Domain_polar_nucleus* polar_nuc = dynamic_cast<const Domain_polar_nucleus*>(&so) ;
       const Domain_polar_shell* polar_shell = dynamic_cast<const Domain_polar_shell*>(&so) ;
       const Domain_polar_compact* polar_compact = dynamic_cast<const Domain_polar_compact*>(&so) ;  
       const Domain_spheric_periodic_nucleus* period_nuc = dynamic_cast<const Domain_spheric_periodic_nucleus*>(&so) ;
       const Domain_spheric_periodic_shell* period_shell = dynamic_cast<const Domain_spheric_periodic_shell*>(&so) ;
       const Domain_spheric_periodic_compact* period_compact = dynamic_cast<const Domain_spheric_periodic_compact*>(&so) ;  
       const Domain_oned_ori* oned_ori = dynamic_cast<const Domain_oned_ori*>(&so) ;
       const Domain_oned_qcq* oned_qcq = dynamic_cast<const Domain_oned_qcq*>(&so) ;
       const Domain_oned_inf* oned_inf = dynamic_cast<const Domain_oned_inf*>(&so) ;  
       const Domain_spheric_time_nucleus* time_nuc = dynamic_cast<const Domain_spheric_time_nucleus*>(&so) ;
       const Domain_spheric_time_shell* time_shell = dynamic_cast<const Domain_spheric_time_shell*>(&so) ;
       const Domain_spheric_time_compact* time_compact = dynamic_cast<const Domain_spheric_time_compact*>(&so) ;
       const Domain_shell_outer_adapted* shell_outer_adapted = dynamic_cast<const Domain_shell_outer_adapted*>(&so) ;
       const Domain_shell_inner_adapted* shell_inner_adapted = dynamic_cast<const Domain_shell_inner_adapted*>(&so) ;
       const Domain_polar_shell_outer_adapted* polar_outer_adapted = dynamic_cast<const Domain_polar_shell_outer_adapted*>(&so) ;
       const Domain_polar_shell_inner_adapted* polar_inner_adapted = dynamic_cast<const Domain_polar_shell_inner_adapted*>(&so) ;
       const Domain_shell_inner_homothetic* shell_inner_homothetic = dynamic_cast<const Domain_shell_inner_homothetic*>(&so) ;
       const Domain_shell_outer_homothetic* shell_outer_homothetic = dynamic_cast<const Domain_shell_outer_homothetic*>(&so) ;
       const Domain_nucleus_symphi* nucsym = dynamic_cast<const Domain_nucleus_symphi*>(&so) ;
       const Domain_shell_symphi* shellsym = dynamic_cast<const Domain_shell_symphi*>(&so) ;
       const Domain_compact_symphi* compsym = dynamic_cast<const Domain_compact_symphi*>(&so) ;
       const Domain_polar_periodic_nucleus* polarperiodicnuc = dynamic_cast<const Domain_polar_periodic_nucleus*>(&so) ;
       const Domain_polar_periodic_shell* polarperiodicshell = dynamic_cast<const Domain_polar_periodic_shell*>(&so) ;
       const Domain_fourD_periodic_nucleus* fourDperiodicnuc = dynamic_cast<const Domain_fourD_periodic_nucleus*>(&so) ;
       const Domain_fourD_periodic_shell* fourDperiodicshell = dynamic_cast<const Domain_fourD_periodic_shell*>(&so) ;

       if (nuc != nullptr)
           o << *nuc << endl ;
       if ((shell != nullptr) && (shell_log==nullptr) && (shell_surr==nullptr))
           o << *shell << endl ;   
       if (shell_log != nullptr) 
           o << *shell_log << endl ; 
       if (shell_surr != nullptr) 
           o << *shell_surr << endl ;
       if (compact !=nullptr) 
           o << *compact << endl ;
       if (bi_rect !=nullptr)
           o << *bi_rect << endl ;
       if (bi_chi !=nullptr)
           o << *bi_chi << endl ;
       if (bi_eta !=nullptr)
           o << *bi_eta << endl ;
	if (critic_inner !=nullptr)
           o << *critic_inner << endl ;
	if (critic_outer !=nullptr)
           o << *critic_outer << endl ;
	if (polar_nuc != nullptr)
           o << *polar_nuc << endl ;
       if (polar_shell != nullptr) 
           o << *polar_shell << endl ;
       if (polar_compact !=nullptr) 
           o << *polar_compact << endl ;
	  if (period_nuc != nullptr)
           o << *period_nuc << endl ;
       if (period_shell != nullptr) 
           o << *period_shell << endl ;
       if (period_compact !=nullptr) 
           o << *period_compact << endl ;
       if (oned_ori != nullptr)
           o << *oned_ori << endl ;
       if (oned_qcq != nullptr) 
           o << *oned_qcq << endl ;
       if (oned_inf !=nullptr) 
           o << *oned_inf << endl ;  
       if (time_nuc != nullptr)
           o << *time_nuc << endl ;
       if (time_shell != nullptr) 
           o << *time_shell << endl ; 
       if (time_compact != nullptr) 
           o << *time_compact << endl ;	
      if (shell_outer_homothetic!=nullptr) 
	o << *shell_outer_homothetic << endl ;
      else {
       if (shell_outer_adapted !=nullptr)
	  o << *shell_outer_adapted << endl ;
      }
	if (shell_inner_homothetic !=nullptr)
	  o << *shell_inner_homothetic << endl ;
	else {
       if (shell_inner_adapted !=nullptr)
	  o << *shell_inner_adapted << endl ;
	}
	
       if (polar_outer_adapted !=nullptr)
	  o << *polar_outer_adapted << endl ;
       if (polar_inner_adapted !=nullptr)
	  o << *polar_inner_adapted << endl ;
      if (nucsym != nullptr)
           o << *nucsym << endl ;
      if (shellsym != nullptr)
           o << *shellsym << endl ;
       if (compsym != nullptr)
           o << *compsym << endl ;
       if (polarperiodicnuc != nullptr)
		o << *polarperiodicnuc << endl ; 
	if (polarperiodicshell != nullptr)
		o << *polarperiodicshell << endl ; 
	if (fourDperiodicnuc != nullptr)
		o << *fourDperiodicnuc << endl ; 
	if (fourDperiodicshell != nullptr)
		o << *fourDperiodicshell << endl ;
       return o ;   
}

// Term eq version of operators
Term_eq Domain::der_normal_term_eq (const Term_eq& so, int bound) const {
  return do_comp_by_comp_with_int (so, bound, &Domain::der_normal) ;
}

Term_eq Domain::dr_term_eq (const Term_eq& so) const {
  return do_comp_by_comp (so, &Domain::der_r) ;
}

Term_eq Domain::dtime_term_eq (const Term_eq& so) const {
  return do_comp_by_comp (so, &Domain::dtime) ;
}

Term_eq Domain::ddtime_term_eq (const Term_eq& so) const {
  return do_comp_by_comp (so, &Domain::ddtime) ;
}

Term_eq Domain::lap_term_eq (const Term_eq& so, int mm) const {
  return do_comp_by_comp_with_int (so, mm, &Domain::laplacian) ;
}

Term_eq Domain::lap2_term_eq (const Term_eq& so, int mm) const {
  return do_comp_by_comp_with_int (so, mm, &Domain::laplacian2) ;
}

Term_eq Domain::mult_r_term_eq (const Term_eq& so) const {
  return do_comp_by_comp (so, &Domain::mult_r) ;
}

Term_eq Domain::div_r_term_eq (const Term_eq& so) const {
  return do_comp_by_comp (so, &Domain::div_r) ;
}

Term_eq Domain::grad_term_eq (const Term_eq& so) const {
 
    Term_eq res (partial_cart(so)) ;
    if (so.val_t->is_m_quant_affected()) {
	  res.val_t->set_parameters().set_m_quant() = so.val_t->get_parameters().get_m_quant() ;
	}
	
    if (so.der_t!=nullptr) {
        if (so.der_t->is_m_quant_affected()) {
            res.der_t->set_parameters().set_m_quant() = so.der_t->get_parameters().get_m_quant() ;
        }
    }
    return res ;
}
Term_eq Domain::integ_term_eq (const Term_eq& so, int bound) const {

	// Check it is a tensor
	if (so.type_data != TERM_T) {
		cerr << "integ_term_eq only defined with respect for a tensor" << endl ;
		abort() ;
	}

	if (so.val_t->get_n_comp() != 1) {
		cerr << "integ_term_eq only defined with respect to a scalar" << endl ;
		abort() ;
	}
	
	assert (so.dom==num_dom) ;
	

	// The value
	Array<int> ind (so.val_t->indices(0)) ;
	Val_domain value ((*so.val_t)(ind)(so.dom)) ;
	double resval ;
	if (value.check_if_zero()) 
		resval= 0. ;
	else 
		resval = integ(value, bound) ;

	if (so.der_t!=nullptr) {
		Val_domain value_der ((*so.der_t)(ind)(so.dom)) ;
		double resder ;
		if (value_der.check_if_zero()) 
			resder = 0. ;
		else 
			resder = integ(value_der, bound) ;
		Term_eq res (so.dom, resval, resder) ;
		return res ;
	}
	else {
		Term_eq res (so.dom, resval) ;
		return res ;
	}
}

Term_eq Domain::do_comp_by_comp_with_int (const Term_eq& target, int val, Val_domain (Domain::*pfunc) (const Val_domain&, int) const ) const {

	// The value	
	Tensor resval (*target.val_t, false) ;

	for (int i=0 ; i<target.val_t->get_n_comp() ; i++) {
		Array<int> ind (target.val_t->indices(i)) ;
		Val_domain value ((*target.val_t)(ind)(num_dom)) ;
		if (value.check_if_zero())
			resval.set(ind).set_domain(num_dom).set_zero() ;
		else {
			Val_domain auxi ((this->*pfunc)(value, val)) ;
			resval.set(ind).set_domain(num_dom) = auxi ;
		}
		if (target.get_val_t().get_parameters())
		  resval.set_parameters() = target.get_val_t().get_parameters() ;

		// affecte parameters 
		if (target.get_val_t().is_name_affected()) {
			resval.set_name_affected() ;
			for (int cmp=0 ; cmp<resval.get_valence() ; cmp++) 
				resval.name_indice[cmp] = target.get_val_t().name_indice[cmp] ;
		}
	}

	if (target.der_t!=nullptr) {
		Tensor resder (*target.der_t, false) ;
		for (int i=0 ; i<target.der_t->get_n_comp() ; i++) {
			Array<int> ind (target.der_t->indices(i)) ;
			Val_domain value ((*target.der_t)(ind)(num_dom)) ;
			if (value.check_if_zero())
				resder.set(ind).set_domain(num_dom).set_zero() ;
			else {
				Val_domain auxi ((this->*pfunc)(value, val)) ;
				resder.set(ind).set_domain(num_dom) = auxi ;
			}
			}
		if (target.get_der_t().get_parameters())
		  resder.set_parameters() = target.get_der_t().get_parameters() ;
		// affecte parameters 
		if (target.get_der_t().is_name_affected()) {
			resder.set_name_affected() ;
			for (int cmp=0 ; cmp<resder.get_valence() ; cmp++) 
				resder.name_indice[cmp] = target.get_der_t().name_indice[cmp] ;
		}

		Term_eq res (num_dom, resval, resder) ;
		return res ;
	}
	else {
		Term_eq res (num_dom, resval) ;
		return res ;
	}
}

Term_eq Domain::do_comp_by_comp (const Term_eq& target, Val_domain (Domain::*pfunc) (const Val_domain&) const ) const {

	// The value	
	Tensor resval (*target.val_t, false) ;

	for (int i=0 ; i<target.val_t->get_n_comp() ; i++) {
		Array<int> ind (target.val_t->indices(i)) ;
		Val_domain value ((*target.val_t)(ind)(num_dom)) ;
		if (value.check_if_zero())
			resval.set(ind).set_domain(num_dom).set_zero() ;
		else {
			Val_domain auxi ((this->*pfunc)(value)) ;
			resval.set(ind).set_domain(num_dom) = auxi ;
		}
		if (target.get_val_t().get_parameters())
		  resval.set_parameters() = target.get_val_t().get_parameters() ;

		// affecte parameters 
		if (target.get_val_t().is_name_affected()) {
			resval.set_name_affected() ;
			for (int cmp=0 ; cmp<resval.get_valence() ; cmp++) 
				resval.name_indice[cmp] = target.get_val_t().name_indice[cmp] ;
		}
	}

	if (target.der_t!=nullptr) {
		Tensor resder (*target.der_t, false) ;
		for (int i=0 ; i<target.der_t->get_n_comp() ; i++) {
			Array<int> ind (target.der_t->indices(i)) ;
			Val_domain value ((*target.der_t)(ind)(num_dom)) ;
			if (value.check_if_zero())
				resder.set(ind).set_domain(num_dom).set_zero() ;
			else {
				Val_domain auxi ((this->*pfunc)(value)) ;
				resder.set(ind).set_domain(num_dom) = auxi ;
			}
			}
		if (target.get_der_t().get_parameters())
		  resder.set_parameters() = target.get_der_t().get_parameters() ;

		// affecte parameters 
		if (target.get_der_t().is_name_affected()) {
			resder.set_name_affected() ;
			for (int cmp=0 ; cmp<resder.get_valence() ; cmp++) 
				resder.name_indice[cmp] = target.get_der_t().name_indice[cmp] ;
		}
		Term_eq res (num_dom, resval, resder) ;
		return res ;
	}
	else {
		Term_eq res (num_dom, resval) ;
		return res ;
	}
}

Term_eq Domain::partial_cart (const Term_eq& so) const {
	int dom = so.get_dom() ;
	int val_res = so.val_t->get_valence() + 1 ;
	Array<int> type_ind (val_res) ;
	type_ind.set(0) = COV ;
	for (int i=1 ; i<val_res ; i++)
		type_ind.set(i) = so.val_t->get_index_type(i-1) ;

	Base_tensor basis (so.get_val_t().get_space()) ;
	basis.set_basis(dom) = CARTESIAN_BASIS ;
	
	// Tensor for val
	Tensor auxi_val (so.get_val_t().get_space(), val_res, type_ind, basis) ;
	
	//Loop on the components :
	Index pos_auxi(auxi_val) ;
	Index pos_so (*so.val_t) ;
	do {
		for (int i=0 ; i<val_res-1 ; i++)
			pos_so.set(i) = pos_auxi(i+1) ;
		auxi_val.set(pos_auxi).set_domain(dom) = (*so.val_t)(pos_so)(dom).der_abs(pos_auxi(0)+1) ;
	}
	while (pos_auxi.inc()) ;

	if (so.der_t==nullptr) {
		return Term_eq (dom, auxi_val) ;
	}
	else {
		// Need to compute the derivative :
		// Tensor for der
		Tensor auxi_der (so.get_val_t().get_space(), val_res, type_ind, basis) ;
		
		//Loop on the components :
		Index pos_auxi_der(auxi_der) ;
		do {
			for (int i=0 ; i<val_res-1 ; i++)
				pos_so.set(i) = pos_auxi_der(i+1) ;
			auxi_der.set(pos_auxi_der).set_domain(dom) = (*so.der_t)(pos_so)(dom).der_abs(pos_auxi_der(0)+1) ;
		}
		while (pos_auxi_der.inc()) ;

		return Term_eq  (dom, auxi_val, auxi_der) ;
	      }
}

Term_eq Domain::partial_spher (const Term_eq& so) const {
  
  
  
	int dom = so.get_dom() ;
	int val_res = so.val_t->get_valence() + 1 ;
	Array<int> type_ind (val_res) ;
	type_ind.set(0) = COV ;
	for (int i=1 ; i<val_res ; i++)
		type_ind.set(i) = so.val_t->get_index_type(i-1) ;

	Base_tensor basis (so.get_val_t().get_space()) ;
	basis.set_basis(dom) = SPHERICAL_BASIS ;
	
	// Tensor for val
	Tensor auxi_val (so.get_val_t().get_space(), val_res, type_ind, basis) ;
	
	//Loop on the components :
	Index pos_auxi(auxi_val) ;
	Index pos_so (*so.val_t) ;
	do {
		for (int i=0 ; i<val_res-1 ; i++)
			pos_so.set(i) = pos_auxi(i+1) ;
		switch (pos_auxi(0)) {
		    case 0 : 
		      // d/dr :
		      auxi_val.set(pos_auxi).set_domain(dom) = (*so.val_t)(pos_so)(dom).der_r() ;
		      break ;
		    case 1 :
		      // 1/r dtheta
		      auxi_val.set(pos_auxi).set_domain(dom) = (*so.val_t)(pos_so)(dom).der_var(2).div_r() ;
		      break ;
		    case 2 :
		      // 1/r sint d/dphi
		      auxi_val.set(pos_auxi).set_domain(dom) = (*so.val_t)(pos_so)(dom).der_var(3).div_r().div_sin_theta() ;
		      break ;
		    default :
		      cerr << "Bad indice in Domain::derive_partial_spher" << endl  ;
		      abort() ;
		}
	}
	while (pos_auxi.inc()) ;

	if (so.der_t==nullptr) {
		return Term_eq (dom, auxi_val) ;
	}
	else {
		// Need to compute the derivative :
		// Tensor for der
		Tensor auxi_der (so.get_val_t().get_space(), val_res, type_ind, basis) ;
		
		//Loop on the components :
		Index pos_auxi_der(auxi_der) ;
		do {
			for (int i=0 ; i<val_res-1 ; i++)
				pos_so.set(i) = pos_auxi_der(i+1) ;
			switch (pos_auxi_der(0)) {
			  case 0 : 
			  // d/dr :
			  auxi_der.set(pos_auxi_der).set_domain(dom) = (*so.der_t)(pos_so)(dom).der_r() ;
			  break ;
			case 1 :
			  // 1/r dtheta
			  auxi_der.set(pos_auxi_der).set_domain(dom) = (*so.der_t)(pos_so)(dom).der_var(2).div_r() ;
			  break ;
			case 2 :
			  // 1/r sint d/dphi
			  auxi_der.set(pos_auxi_der).set_domain(dom) = (*so.der_t)(pos_so)(dom).der_var(3).div_r().div_sin_theta() ;
			  break ;
			default :
			  cerr << "Bad indice in Domain::derive_partial_spher" << endl  ;
			  abort() ;
		    }
		}
		while (pos_auxi_der.inc()) ;

		return Term_eq  (dom, auxi_val, auxi_der) ;
	      }
}

Term_eq Domain::partial_mtz (const Term_eq& so) const {
  
  
  
	int dom = so.get_dom() ;
	int val_res = so.val_t->get_valence() + 1 ;
	Array<int> type_ind (val_res) ;
	type_ind.set(0) = COV ;
	for (int i=1 ; i<val_res ; i++)
		type_ind.set(i) = so.val_t->get_index_type(i-1) ;

	Base_tensor basis (so.get_val_t().get_space()) ;
	basis.set_basis(dom) = MTZ_BASIS ;
	
	// Tensor for val
	Tensor auxi_val (so.get_val_t().get_space(), val_res, type_ind, basis) ;
	
	//Loop on the components :
	Index pos_auxi(auxi_val) ;
	Index pos_so (*so.val_t) ;
	do {
		for (int i=0 ; i<val_res-1 ; i++)
			pos_so.set(i) = pos_auxi(i+1) ;
		switch (pos_auxi(0)) {
		    case 0 : 
		      // d/dr :
		      auxi_val.set(pos_auxi).set_domain(dom) = (*so.val_t)(pos_so)(dom).der_r() ;
		      break ;
		    case 1 :
		      // cost/r dtheta
		      auxi_val.set(pos_auxi).set_domain(dom) = (*so.val_t)(pos_so)(dom).der_var(2).div_r().mult_cos_theta() ;
		      break ;
		    case 2 :
		      // cost/r sint d/dphi
		      auxi_val.set(pos_auxi).set_domain(dom) = (*so.val_t)(pos_so)(dom).der_var(3).div_r().div_sin_theta().mult_cos_theta() ;
		      break ;
		    default :
		      cerr << "Bad indice in Domain::derive_partial_mtz" << endl  ;
		      abort() ;
		}
	}
	while (pos_auxi.inc()) ;

	if (so.der_t==nullptr) {
		return Term_eq (dom, auxi_val) ;
	}
	else {
		// Need to compute the derivative :
		// Tensor for der
		Tensor auxi_der (so.get_val_t().get_space(), val_res, type_ind, basis) ;
		
		//Loop on the components :
		Index pos_auxi_der(auxi_der) ;
		do {
			for (int i=0 ; i<val_res-1 ; i++)
				pos_so.set(i) = pos_auxi_der(i+1) ;
			switch (pos_auxi_der(0)) {
			  case 0 : 
			  // d/dr :
			  auxi_der.set(pos_auxi_der).set_domain(dom) = (*so.der_t)(pos_so)(dom).der_r() ;
			  break ;
			case 1 :
			  // 1/r dtheta
			  auxi_der.set(pos_auxi_der).set_domain(dom) = (*so.der_t)(pos_so)(dom).der_var(2).div_r().mult_cos_theta() ;
			  break ;
			case 2 :
			  // 1/r sint d/dphi
			  auxi_der.set(pos_auxi_der).set_domain(dom) = (*so.der_t)(pos_so)(dom).der_var(3).div_r().div_sin_theta().mult_cos_theta() ;
			  break ;
			default :
			  cerr << "Bad indice in Domain::derive_partial_mtz" << endl  ;
			  abort() ;
		    }
		}
		while (pos_auxi_der.inc()) ;

		return Term_eq  (dom, auxi_val, auxi_der) ;
	      }
}


Term_eq Domain::connection_spher (const Term_eq& so) const {
  
  int dom = so.get_dom() ;
  assert (dom == num_dom) ;
  
  int valence = so.val_t->get_valence() ;
  int val_res = so.val_t->get_valence() + 1 ;

  Array<int> type_ind (val_res) ;
  type_ind.set(0) = COV ;
  for (int i=1 ; i<val_res ; i++)
      type_ind.set(i) = so.val_t->get_index_type(i-1) ;
  
  Base_tensor basis (so.get_val_t().get_space()) ;
  basis.set_basis(dom) = SPHERICAL_BASIS ;
	
  Tensor auxi_val (so.get_val_t().get_space(), val_res, type_ind, basis, 3) ;
  for (int cmp=0 ; cmp<auxi_val.get_n_comp() ; cmp ++)
    auxi_val.set(auxi_val.indices(cmp)) = 0  ;
  
  for (int ind_sum=0 ; ind_sum<valence ; ind_sum++) {

	  //Loop on the components :
	  Index pos_auxi(auxi_val) ;
	  Index pos_so (*so.val_t) ;

	  do {
		  for (int i=0 ; i<valence ; i++)
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
			      auxi_val.set(pos_auxi).set_domain(dom) -= (*so.val_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 1 :
			      // Dtheta S_theta
			      pos_so.set(ind_sum) = 0 ;
			      auxi_val.set(pos_auxi).set_domain(dom) += (*so.val_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 2 :
			      //Dtheta S_phi 
			      break ;
			    default :
			      cerr << "Bad indice in Domain::connection_spher" << endl  ;
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
			      auxi_val.set(pos_auxi).set_domain(dom) -= (*so.val_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 1 :
			      // Dphi S_theta
			      pos_so.set(ind_sum) = 2 ;
			      auxi_val.set(pos_auxi).set_domain(dom) -= (*so.val_t)(pos_so)(dom).div_r().mult_cos_theta().div_sin_theta() ;
			      break ;
			    case 2 :
			      //Dphi S_phi
			      pos_so.set(ind_sum) = 0 ;
			      auxi_val.set(pos_auxi).set_domain(dom) += (*so.val_t)(pos_so)(dom).div_r() ;
			      pos_so.set(ind_sum) = 1 ;
			      auxi_val.set(pos_auxi).set_domain(dom) += (*so.val_t)(pos_so)(dom).div_r().mult_cos_theta().div_sin_theta() ;
			      break ;
			    default :
			      cerr << "Bad indice in Domain::connection_spher" << endl  ;
			      abort() ;
			}
			break ;
		      default :
			cerr << "Bad indice in Domain::connection_spher" << endl  ;
			abort() ;
		  }
	  }
	  while (pos_auxi.inc()) ;
	}
	
	
	if (so.der_t==nullptr) {
		// No need for derivative :
		return Term_eq (dom, auxi_val) ;
	}
	else {
	
	  // Need to compute the derivative :
	  // Tensor for der
	  Tensor auxi_der (so.get_val_t().get_space(), val_res, type_ind, basis, 3) ;
	  for (int cmp=0 ; cmp<auxi_der.get_n_comp() ; cmp ++)
	    auxi_der.set(auxi_der.indices(cmp)) = 0  ;
 
		// Loop indice summation on connection symbols 
	for (int ind_sum=0 ; ind_sum<valence ; ind_sum++) {

	  //Loop on the components :
	  Index pos_auxi_der(auxi_der) ;
	  Index pos_so (*so.der_t) ;

	  do {
		  for (int i=0 ; i<valence ; i++)
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
			      auxi_der.set(pos_auxi_der).set_domain(dom) -= (*so.der_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 1 :
			      // Dtheta S_theta
			      pos_so.set(ind_sum) = 0 ;
			      auxi_der.set(pos_auxi_der).set_domain(dom) += (*so.der_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 2 :
			      //Dtheta S_phi 
			      break ;
			    default :
			      cerr << "Bad indice in Domain::connection_spher" << endl  ;
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
			      auxi_der.set(pos_auxi_der).set_domain(dom) -= (*so.der_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 1 :
			      // Dphi S_theta
			      pos_so.set(ind_sum) = 2 ;
			      auxi_der.set(pos_auxi_der).set_domain(dom) -= (*so.der_t)(pos_so)(dom).div_r().mult_cos_theta().div_sin_theta() ;
			      break ;
			    case 2 :
			      //Dphi S_phi
			      pos_so.set(ind_sum) = 0 ;
			      auxi_der.set(pos_auxi_der).set_domain(dom) += (*so.der_t)(pos_so)(dom).div_r() ;
			      pos_so.set(ind_sum) = 1 ;
			      auxi_der.set(pos_auxi_der).set_domain(dom) += (*so.der_t)(pos_so)(dom).div_r().mult_cos_theta().div_sin_theta() ;
			      break ;
			    default :
			      cerr << "Bad indice in Domain::connection_spher" << endl  ;
			      abort() ;
			}
			break ;
		      default :
			cerr << "Bad indice in Domain::connection_spher" << endl  ;
			abort() ;
		  }
	  }
	  while (pos_auxi_der.inc()) ;
	}

	return Term_eq (dom, auxi_val, auxi_der) ;
	}
}

Term_eq Domain::connection_mtz (const Term_eq& so) const {
  
  int dom = so.get_dom() ;
  assert (dom == num_dom) ;
  
  int valence = so.val_t->get_valence() ;
  int val_res = so.val_t->get_valence() + 1 ;

  Array<int> type_ind (val_res) ;
  type_ind.set(0) = COV ;
  for (int i=1 ; i<val_res ; i++)
      type_ind.set(i) = so.val_t->get_index_type(i-1) ;
  
  Base_tensor basis (so.get_val_t().get_space()) ;
  basis.set_basis(dom) = MTZ_BASIS ;
	
  Tensor auxi_val (so.get_val_t().get_space(), val_res, type_ind, basis, 3) ;
  for (int cmp=0 ; cmp<auxi_val.get_n_comp() ; cmp ++)
    auxi_val.set(auxi_val.indices(cmp)) = 0  ;
  
  for (int ind_sum=0 ; ind_sum<valence ; ind_sum++) {

	  //Loop on the components :
	  Index pos_auxi(auxi_val) ;
	  Index pos_so (*so.val_t) ;

	  do {
		  for (int i=0 ; i<valence ; i++)
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
			      auxi_val.set(pos_auxi).set_domain(dom) -= (*so.val_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 1 :
			      // Dtheta S_theta
			      pos_so.set(ind_sum) = 0 ;
			      auxi_val.set(pos_auxi).set_domain(dom) += (*so.val_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 2 :
			      //Dtheta S_phi 
			      break ;
			    default :
			      cerr << "Bad indice in Domain::connection_mtz" << endl  ;
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
			      auxi_val.set(pos_auxi).set_domain(dom) -= (*so.val_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 1 :
			      // Dphi S_theta
			      pos_so.set(ind_sum) = 2 ;
			      auxi_val.set(pos_auxi).set_domain(dom) -= (*so.val_t)(pos_so)(dom).div_r().div_sin_theta() ;
			      break ;
			    case 2 :
			      //Dphi S_phi
			      pos_so.set(ind_sum) = 0 ;
			      auxi_val.set(pos_auxi).set_domain(dom) += (*so.val_t)(pos_so)(dom).div_r() ;
			      pos_so.set(ind_sum) = 1 ;
			      auxi_val.set(pos_auxi).set_domain(dom) += (*so.val_t)(pos_so)(dom).div_r().div_sin_theta() ;
			      break ;
			    default :
			      cerr << "Bad indice in Domain::connection_mtz" << endl  ;
			      abort() ;
			}
			break ;
		      default :
			cerr << "Bad indice in Domain::connection_mtz" << endl  ;
			abort() ;
		  }
	  }
	  while (pos_auxi.inc()) ;
	}
	
	
	if (so.der_t==nullptr) {
		// No need for derivative :
		return Term_eq (dom, auxi_val) ;
	}
	else {
	
	  // Need to compute the derivative :
	  // Tensor for der
	  Tensor auxi_der (so.get_val_t().get_space(), val_res, type_ind, basis, 3) ;
	  for (int cmp=0 ; cmp<auxi_der.get_n_comp() ; cmp ++)
	    auxi_der.set(auxi_der.indices(cmp)) = 0  ;
 
		// Loop indice summation on connection symbols 
	for (int ind_sum=0 ; ind_sum<valence ; ind_sum++) {

	  //Loop on the components :
	  Index pos_auxi_der(auxi_der) ;
	  Index pos_so (*so.der_t) ;

	  do {
		  for (int i=0 ; i<valence ; i++)
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
			      auxi_der.set(pos_auxi_der).set_domain(dom) -= (*so.der_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 1 :
			      // Dtheta S_theta
			      pos_so.set(ind_sum) = 0 ;
			      auxi_der.set(pos_auxi_der).set_domain(dom) += (*so.der_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 2 :
			      //Dtheta S_phi 
			      break ;
			    default :
			      cerr << "Bad indice in Domain::connection_mtz" << endl  ;
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
			      auxi_der.set(pos_auxi_der).set_domain(dom) -= (*so.der_t)(pos_so)(dom).div_r() ;
			      break ;
			    case 1 :
			      // Dphi S_theta
			      pos_so.set(ind_sum) = 2 ;
			      auxi_der.set(pos_auxi_der).set_domain(dom) -= (*so.der_t)(pos_so)(dom).div_r().div_sin_theta() ;
			      break ;
			    case 2 :
			      //Dphi S_phi
			      pos_so.set(ind_sum) = 0 ;
			      auxi_der.set(pos_auxi_der).set_domain(dom) += (*so.der_t)(pos_so)(dom).div_r() ;
			      pos_so.set(ind_sum) = 1 ;
			      auxi_der.set(pos_auxi_der).set_domain(dom) += (*so.der_t)(pos_so)(dom).div_r().div_sin_theta() ;
			      break ;
			    default :
			      cerr << "Bad indice in Domain::connection_mtz" << endl  ;
			      abort() ;
			}
			break ;
		      default :
			cerr << "Bad indice in Domain::connection_mtz" << endl  ;
			abort() ;
		  }
	  }
	  while (pos_auxi_der.inc()) ;
	}

	return Term_eq (dom, auxi_val, auxi_der) ;
	}
}


Term_eq Domain::integ_volume_term_eq (const Term_eq& target) const {
  
	int dom = target.get_dom() ;
  
	// Check it is a tensor
	if (target.type_data != TERM_T) {
		cerr << "Ope_int_volume only defined with respect for a tensor" << endl ;
		abort() ;
	}

	if (target.val_t->get_n_comp() != 1) {
		cerr << "Ope_int_volume only defined with respect to a scalar" << endl ;
		abort() ;
	}

	// The value
	Array<int> ind (target.val_t->indices(0)) ;
	Val_domain value ((*target.val_t)(ind)(dom)) ;
	double resval ;
	if (value.check_if_zero()) 
		resval= 0. ;
	else 
		resval = value.get_domain()->integ_volume(value) ;

	if (target.der_t!=nullptr) {
		Val_domain value_der ((*target.der_t)(ind)(dom)) ;
		double resder ;
		if (value_der.check_if_zero()) 
			resder = 0. ;
		else 
			resder = value_der.get_domain()->integ_volume(value_der) ;
		Term_eq res (dom, resval, resder) ;
		return res ;
	}
	else {
		Term_eq res (dom, resval) ;
		return res ;
	}
}

// List of all the functions that must be implemented for the concerned types of domain
bool Domain::is_in(const Point&, double) const {
	cerr << "is_in not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

const Point Domain::absol_to_num(const Point&) const {
	cerr << "Absol_to_num not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

const Point Domain::absol_to_num_bound(const Point&, int) const {
	cerr << "Absol_to_num_bound not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::do_der_abs_from_der_var(const Val_domain *const *const der_var, Val_domain **const der_abs) const {
	cerr << "do_der_abs_from_der_var not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}
   
Base_spectral Domain::mult (const Base_spectral&, const Base_spectral&) const {
	cerr << "mult not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::do_coloc () {
	cerr << "do_coloc not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_cheb_base(Base_spectral&) const {
	cerr << "Symetric Chebyshev spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}
void Domain::set_legendre_base(Base_spectral&) const {
	cerr << "Symetric Legendre spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_cheb_r_base(Base_spectral&) const {
	cerr << "Chebyshev spectral base for r not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_legendre_r_base(Base_spectral&) const {
	cerr << "Legendre spectral base for r not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_anti_cheb_base(Base_spectral&) const {
	cerr << "Anti symetric Chebyshev spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_anti_legendre_base(Base_spectral&) const {
	cerr << "Anti symetric Legendre spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}
void Domain::set_cheb_base_with_m(Base_spectral&, int) const {
	cerr << "Symetric Chebyshev spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_legendre_base_with_m(Base_spectral&, int) const {
	cerr << "Symetric Legendre spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_anti_cheb_base_with_m(Base_spectral&, int) const {
	cerr << "Anti symetric Chebyshev spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_anti_legendre_base_with_m(Base_spectral&, int) const {
	cerr << "Anti symetric Legendre spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_cheb_xodd_base(Base_spectral&) const {
	cerr << "Chebyshev with xodd spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_legendre_xodd_base(Base_spectral&) const {
	cerr << "Legendre with xodd spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_cheb_todd_base(Base_spectral&) const {
	cerr << "Chebyshev with todd spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_legendre_todd_base(Base_spectral&) const {
	cerr << "Legendre with todd spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_cheb_xodd_todd_base(Base_spectral&) const {
	cerr << "Chebyshev with X and T odd spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_legendre_xodd_todd_base(Base_spectral&) const {
	cerr << "Odd Legendre with X and T odd spectral base not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_cheb_base_r_spher (Base_spectral&) const {
    cerr << "Cheb base r not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_cheb_base_t_spher (Base_spectral&) const {
    cerr << "Cheb base t not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_cheb_base_p_spher (Base_spectral&) const {
    cerr << "Cheb base p not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}
void Domain::set_cheb_base_r_mtz (Base_spectral&) const {
    cerr << "Cheb base r MTZ not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_cheb_base_t_mtz  (Base_spectral&) const {
    cerr << "Cheb base t MTZ not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_cheb_base_p_mtz (Base_spectral&) const {
    cerr << "Cheb base p MTZ not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}
void Domain::set_cheb_base_rt_spher (Base_spectral&) const {
    cerr << "Cheb base rt not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_cheb_base_rp_spher (Base_spectral&) const {
    cerr << "Cheb base rp not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_cheb_base_tp_spher (Base_spectral&) const {
    cerr << "Cheb base tp not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}
void Domain::set_legendre_base_r_spher (Base_spectral&) const {
    cerr << "Legendre base r not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_legendre_base_t_spher (Base_spectral&) const {
    cerr << "Legendre base t not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_legendre_base_p_spher (Base_spectral&) const {
    cerr << "Legendre base p not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_legendre_base_r_mtz (Base_spectral&) const {
    cerr << "Legendre base r MTZ not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_legendre_base_t_mtz (Base_spectral&) const {
    cerr << "Legendre base t MTZ not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_legendre_base_p_mtz (Base_spectral&) const {
    cerr << "Legendre base p MTZ not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}
void Domain::set_cheb_base_odd (Base_spectral&) const {
    cerr << "Cheb base with odd cosines not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}

void Domain::set_legendre_base_odd (Base_spectral&) const {
    cerr << "Legendre base with odd cosines not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}
void Domain::set_cheb_base_xy_cart(Base_spectral&) const {
    cerr << "Cheb base xy not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}


void Domain::set_cheb_base_xz_cart(Base_spectral&) const {
    cerr << "Cheb base xz not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}


void Domain::set_cheb_base_yz_cart(Base_spectral&) const {
    cerr << "Cheb base yz not implemented for " << endl ;
    cerr << *this << endl ;
    abort() ;
}
void Domain::set_cheb_base_x_cart(Base_spectral& ba) const {
    // Default version
    set_cheb_base(ba) ;
}


void Domain::set_cheb_base_y_cart(Base_spectral& ba) const {
     // Default version
    set_cheb_base(ba) ;
}


void Domain::set_cheb_base_z_cart(Base_spectral& ba) const {
    // Default version
    set_anti_cheb_base(ba) ;
}


void Domain::set_legendre_base_x_cart(Base_spectral& ba) const {
   // Default version
    set_legendre_base(ba) ;
}


void Domain::set_legendre_base_y_cart(Base_spectral& ba) const {
   // Default version
    set_legendre_base(ba) ;
}


void Domain::set_legendre_base_z_cart(Base_spectral& ba) const {
     // Default version
    set_anti_legendre_base(ba) ;
}

  
Val_domain Domain::mult_cos_theta (const Val_domain&) const {
	cerr << "Multiplication by cos(theta) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::mult_cos_phi (const Val_domain&) const {
	cerr << "Multiplication by cos(phi) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::mult_sin_theta (const Val_domain&) const {
	cerr << "Multiplication by sin(theta) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::mult_sin_phi (const Val_domain&) const {
	cerr << "Multiplication by sin(phi) not implemented for this type of domain" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::div_sin_theta (const Val_domain&) const {
	cerr << "Division by sin(theta) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::div_cos_theta (const Val_domain&) const {
	cerr << "Division by cos(theta) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::div_x (const Val_domain&) const {
	cerr << "Division by x not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::set_val_inf (Val_domain&, double) const {
	cerr << "set_val_inf not defined for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::div_xm1 (const Val_domain&) const {
	cerr << "Division by (x-1) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::div_xp1 (const Val_domain&) const {
	cerr << "Division by (x+1) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::div_1mrsL (const Val_domain&) const {
	cerr << "Division by 1 - r/L not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}


Val_domain Domain::mult_1mrsL (const Val_domain&) const {
	cerr << "Multiplication by 1 - r/L  not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::mult_xm1 (const Val_domain&) const {
	cerr << "Multiplication by (x-1) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::mult_r (const Val_domain&) const {
	cerr << "Multiplication by r not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::mult_x (const Val_domain&) const {
	cerr << "Multiplication by x not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::mult_cos_time (const Val_domain&) const {
	cerr << "Multiplication by cos(time) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::mult_sin_time (const Val_domain&) const {
	cerr << "Multiplication by sin(time) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}
   
Val_domain Domain::div_1mx2 (const Val_domain&) const {
	cerr << "Divison by 1-x^2 not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::srdr (const Val_domain&) const {
	cerr << "srdr not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::div_r (const Val_domain&) const {
	cerr << "Division by r not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::div_sin_chi (const Val_domain&) const {
	cerr << "Division by sin(chi) not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Val_domain Domain::div_chi (const Val_domain&) const {
	cerr << "Division by chi not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Tensor Domain::change_basis_cart_to_spher (int, const Tensor&) const {
	cerr << "change_basis_cart_to_spher not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Tensor Domain::change_basis_spher_to_cart (int, const Tensor&) const {
	cerr << "change_basis_spher_to_cart not implemented for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

double Domain::get_rmax() const {
	cerr << "No rmax for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

double Domain::get_rmin() const {
	cerr << "No rmin for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

Point Domain::get_center () const {
	cerr << "No center for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

const Val_domain & Domain::get_chi() const {
	cerr << "No chi for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

const Val_domain & Domain::get_eta() const {
	cerr << "No eta for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

const Val_domain & Domain::get_X() const {
	cerr << "No X for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

const Val_domain & Domain::get_T() const {
	cerr << "No T for" << endl ;
	cerr << *this << endl  ;
	abort() ;
}

void Domain::find_other_dom (int, int, int&, int&) const {
	cerr << "find_other_dom not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::der_normal (const Val_domain&, int) const {
	cerr << "der_normal not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::der_partial_var (const Val_domain&, int) const {
	cerr << "der_partial_var not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::der_r (const Val_domain&) const {
	cerr << "der_r not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::der_t (const Val_domain&) const {
	cerr << "der_t not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}
Val_domain Domain::der_p (const Val_domain&) const {
	cerr << "der_p not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::der_r_rtwo (const Val_domain&) const {
	cerr << "der_r_rtwo not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::ddr (const Val_domain& so) const {
	return so.der_r().der_r() ;
}

Val_domain Domain::ddp (const Val_domain&) const {
	cerr << "ddp not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::ddt (const Val_domain&) const {
	cerr << "ddt not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::dt (const Val_domain&) const {
	cerr << "dt not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::dtime (const Val_domain&) const {
	cerr << "dtime not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Val_domain Domain::ddtime(const Val_domain&) const {
	cerr << "ddtime not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

double Domain::integ (const Val_domain&, int) const {
	cerr << "integ not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

double Domain::integrale (const Val_domain&) const {
	cerr << "integrale not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

int Domain::nbr_unknowns (const Tensor&, int) const {
	cerr << "nbr_unknowns not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Array<int> Domain::nbr_conditions (const Tensor&, int, int, int, Array<int>**) const {
	cerr << "nbr_conditions not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Array<int> Domain::nbr_conditions_array (const Tensor&, int, const Array<int>&, int, Array<int>**) const {
	cerr << "nbr_conditions (with array) not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Array<int> Domain::nbr_conditions_boundary (const Tensor&, int, int, int, Array<int>**) const {
	cerr << "nbr_conditions_boundary not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Array<int> Domain::nbr_conditions_boundary_array (const Tensor&, int, int, const Array<int>&, int, Array<int>**) const {
	cerr << "nbr_conditions_boundary (array version) not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}
Array<int> Domain::nbr_conditions_boundary_one_side (const Tensor&, int, int, int, Array<int>**) const {
	cerr << "nbr_conditions_boundary_one_side not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int, Array<int>**) const {
	cerr << "export_tau not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::export_tau_array (const Tensor&, int, const Array<int>&, Array<double>&, int&, const Array<int>&, int, Array<int>**) const {
	cerr << "export_tau (with array) not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::export_tau_boundary (const Tensor&, int, int,  Array<double>&, int&, const Array<int>&, int, Array<int>**) const {
	cerr << "export_tau_boundary not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}
void Domain::export_tau_boundary_exception (const Tensor&, int, int,  Array<double>&, int&, const Array<int>&, const Param&, int, const Tensor&, 
															    int, Array<int>**) const {
	cerr << "export_tau_boundary_exception not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::export_tau_boundary_array (const Tensor&, int, int, const Array<int>&, Array<double>&, int&, const Array<int>&, int, Array<int>**) const {
	cerr << "export_tau_boundary (with array) not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::export_tau_boundary_one_side (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int, Array<int>**) const {
	cerr << "export_tau_boundary_one_side not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::affecte_tau (Tensor&, int, const Array<double>&, int&) const {
	cerr << "affecte_tau not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::affecte_tau_one_coef (Tensor&, int, int, int&) const {
	cerr << "affecte_tau_one_coef not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

double Domain::val_boundary (int, const Val_domain&, const Index&) const {
	cerr << "val_boundary not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

int Domain::nbr_points_boundary (int, const Base_spectral&) const {
	cerr << "nbr_points_boundary not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::do_which_points_boundary (int, const Base_spectral&, Index**, int) const {
	cerr << "do_which_points_boundary not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::do_absol () const {
	cerr << "do_absol not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::do_cart () const {
	cerr << "do_cart not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::do_cart_surr () const {
	cerr << "do_cart_surr not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

void Domain::do_radius () const {
	cerr << "do_radius not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

double Domain::multipoles_sym (int, int, int, const Val_domain&, const Array<double>&) const {
	cerr << "multipoles sym not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

double Domain::multipoles_asym (int, int, int, const Val_domain&, const Array<double>&) const {
	cerr << "multipoles asym not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Term_eq Domain::multipoles_sym (int, int, int, const Term_eq&, const Array<double>&) const {
	cerr << "multipoles sym not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Term_eq Domain::multipoles_asym (int, int, int, const Term_eq&, const Array<double>&) const {
	cerr << "multipoles asym not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Term_eq Domain::radial_part_sym (const Space&, int, int, const Term_eq&, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param&) const {
	cerr << "radial part sym not implemented for" << endl ;
	cerr << *this << endl ;	
	cerr << "and function " << f << endl ;
	abort() ;
}

Term_eq Domain::radial_part_asym (const Space&, int, int, const Term_eq&, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param&) const {
	cerr << "radial part asym not implemented for" << endl ;
	cerr << *this << endl ;	
	cerr << "and function " << f << endl ;
	abort() ;
}

Term_eq Domain::harmonics_sym (const Term_eq&, const Term_eq&, int, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param&, const Array<double>&) const {
	cerr << "harmonics_sym not implemented for" << endl ;
	cerr << *this << endl ;	
	cerr << "and function " << f << endl ;
	abort() ;
}

Term_eq Domain::harmonics_asym (const Term_eq&, const Term_eq&, int, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param&, const Array<double>&) const {
	cerr << "harmonics_asym not implemented for" << endl ;
	cerr << *this << endl ;
	cerr << "and function " << f << endl ;
	abort() ;
}

Term_eq Domain::der_multipoles_sym (int, int, int, const Term_eq&, const Array<double>&) const {
	cerr << "der multipoles sym not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Term_eq Domain::der_multipoles_asym (int, int, int, const Term_eq&, const Array<double>&) const {
	cerr << "der multipoles asym not implemented for" << endl ;
	cerr << *this << endl ;
	abort() ;
}

Term_eq Domain::der_radial_part_sym (const Space&, int, int, const Term_eq&, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param&) const {
	cerr << "der radial part sym not implemented for" << endl ;
	cerr << *this << endl ;	
	cerr << "and function " << f << endl ;
	abort() ;
}

Term_eq Domain::der_radial_part_asym (const Space&, int, int, const Term_eq&, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param&) const {
	cerr << "der radial part asym not implemented for" << endl ;
	cerr << *this << endl ;	
	cerr << "and function " << f << endl ;
	abort() ;
}

Term_eq Domain::der_harmonics_sym (const Term_eq&, const Term_eq&, int, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param&, const Array<double>&) const {
	cerr << "der_harmonics_sym not implemented for" << endl ;
	cerr << *this << endl ;	
	cerr << "and function " << f << endl ;
	abort() ;
}

Term_eq Domain::der_harmonics_asym (const Term_eq&, const Term_eq&, int, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param&, const Array<double>&) const {
	cerr << "der harmonics_asym not implemented for" << endl ;
	cerr << *this << endl ;
	cerr << "and function " << f << endl ;
	abort() ;
}

const Term_eq* Domain::give_normal (int, int) const {
  cerr << "give_normal not implemented for" << endl ;
  cerr << *this << endl ;
  abort() ;
}

Term_eq Domain::derive_flat_spher (int, char, const Term_eq&, const Metric*) const {
    cerr << "derive_flat_spher not implemented for" << endl ;
    cerr << *this << endl ;
    abort() ;
}

Term_eq Domain::derive_flat_mtz (int, char, const Term_eq&, const Metric*) const {
    cerr << "derive_flat_mtz not implemented for" << endl ;
    cerr << *this << endl ;
    abort() ;
}

Term_eq Domain::derive_flat_cart (int, char, const Term_eq&, const Metric*) const {
    cerr << "derive_flat_cart not implemented for" << endl ;
    cerr << *this << endl ;
    abort() ;
}

int Domain::give_place_var (char*) const {
  return -1 ;
}

Tensor Domain::import (int, int, int, const Array<int>&, Tensor**) const {
  cerr << "import not implemented for" << endl ;
  cerr << *this << endl ;
  abort() ;
}

void Domain::update_term_eq (Term_eq* so) const {
   so->set_der_zero() ;
}

double Domain::integ_volume (const Val_domain&) const {
  cerr << "Integ volume not implemented for " << endl ;
  cerr << *this << endl ;
  abort() ;
}

void Domain::filter (Tensor&, int, double) const {
   cerr << "Filter not implemented for" << endl ;
   cerr << *this << endl ;
   abort() ;
}

Term_eq Domain::div_1mx2_term_eq (const Term_eq& so) const  {
	return do_comp_by_comp (so, &Domain::div_1mx2) ;
}
}
