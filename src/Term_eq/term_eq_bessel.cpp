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

#include "term_eq.hpp"
#include "scalar.hpp"
#include <gsl/gsl_sf_bessel.h>
namespace Kadath {
Term_eq bessel_jl (const Term_eq& so, int ll) {

  int dom = so.get_dom() ;
  
  // Double case 
  if (so.type_data == TERM_D) { 
    double val = gsl_sf_bessel_jl (ll, *so.val_d) ;
    bool doder = (so.der_d == 0x0) ? false : true ;
    if (doder) {
      double der = (ll==0) ? (cos (*so.val_d)/(*so.val_d) - sin(*so.val_d)/(*so.val_d)/(*so.val_d))*(*so.der_d) : 
			     (gsl_sf_bessel_jl (ll-1,(*so.val_d)) - (ll+1)/(*so.val_d)*val)*(*so.der_d)  ;
      return Term_eq (dom, val, der) ;
    }
    else {
      return Term_eq (dom, val) ;
    }
  }
 
   // Scalar case
  if (so.type_data == TERM_T) {
    if (so.val_t->get_valence()!=0) {
      cerr << "Bessel function only defined for scalars so far..." << endl ;
      abort() ;
    }

    Val_domain val (so.val_t->get_space().get_domain(dom)) ;
    val.allocate_conf() ;
    Val_domain der (so.val_t->get_space().get_domain(dom)) ;
    der.allocate_conf() ;
  
    bool doder = (so.der_t == 0x0) ? false : true ;
    Index pos (val.get_domain()->get_nbr_points()) ;
    do  {
      double vv = (*so.val_t)()(dom)(pos) ;
      val.set(pos) = gsl_sf_bessel_jl (ll, vv) ;
    
      if (doder) {
	  double dd = (*so.der_t)()(dom)(pos) ;
	  if (ll==0)
	    der.set(pos) = (cos (vv)/vv - sin(vv)/vv/vv)*dd ;
	  else
	    der.set(pos) = (gsl_sf_bessel_jl (ll-1,vv) - (ll+1)/vv*val(pos))*dd ;
      }
    }
    while (pos.inc()) ;
  
    val.std_base() ;
    if (doder)
      der.std_base() ;
  
    Scalar valres (so.val_t->get_space()) ;
    valres.set_domain(dom) = val ;
    if (doder) {
      Scalar derres (so.val_t->get_space()) ;
      derres.set_domain(dom) = der ;
      return Term_eq (dom, valres, derres) ;
    }
    else  {
      return Term_eq (dom, valres) ;
    }
  }
  
  // Unknown case
  cerr << "Unknown type of Term_eq in bessel_jl" << endl  ;
  abort() ;
}


Term_eq bessel_yl (const Term_eq& so, int ll) {
    
  int dom = so.get_dom() ;
  
  // Double case 
  if (so.type_data == TERM_D) { 
    double val = gsl_sf_bessel_yl (ll, *so.val_d) ;
    bool doder = (so.der_d == 0x0) ? false : true ;
    if (doder) {
      double der = (ll==0) ? (sin (*so.val_d)/(*so.val_d) + cos(*so.val_d)/(*so.val_d)/(*so.val_d))*(*so.der_d) : 
			     (gsl_sf_bessel_yl (ll-1,(*so.val_d)) - (ll+1)/(*so.val_d)*val)*(*so.der_d)  ;
      return Term_eq (dom, val, der) ;
    }
    else {
      return Term_eq (dom, val) ;
    }
  }
  
  // Scalar case
  if (so.type_data == TERM_T) {
    if (so.val_t->get_valence()!=0) {
      cerr << "Bessel function only defined for scalars so far..." << endl ;
      abort() ;
    }

    Val_domain val (so.val_t->get_space().get_domain(dom)) ;
    val.allocate_conf() ;
    Val_domain der (so.val_t->get_space().get_domain(dom)) ;
    der.allocate_conf() ;
  
    bool doder = (so.der_t == 0x0) ? false : true ;
    Index pos (val.get_domain()->get_nbr_points()) ;
    do  {
      double vv = (*so.val_t)()(dom)(pos) ;
      val.set(pos) = gsl_sf_bessel_yl (ll, vv) ;
      if (doder) {
	double dd = (*so.der_t)()(dom)(pos) ;
	if (ll==0)
	    der.set(pos) = (sin (vv)/vv + cos(vv)/vv/vv)*dd ;
	else
	    der.set(pos) = (gsl_sf_bessel_yl (ll-1,vv) - (ll+1)/vv*val(pos))*dd ;
      }
    }
    while (pos.inc()) ;
  
    val.std_base() ;
    if (doder)
      der.std_base() ;
  
    Scalar valres (so.val_t->get_space()) ;
    valres.set_domain(dom) = val ;
    if (doder) {
	Scalar derres (so.val_t->get_space()) ;
	derres.set_domain(dom) = der ;
	return Term_eq (dom, valres, derres) ;
    }
    else  {
	return Term_eq (dom, valres) ;
    }
  }
  
   // Unknown case
  cerr << "Unknown type of Term_eq in bessel_yl" << endl  ;
  abort() ;
}

Term_eq bessel_djl (const Term_eq& so, int ll) {

  int dom = so.get_dom() ;
  
  // Double case 
  if (so.type_data == TERM_D) { 
    double val ;
    
    double jl = gsl_sf_bessel_jl (ll,(*so.val_d)) ;
    double jlm1 = (ll==0) ? 0 : gsl_sf_bessel_jl (ll-1 ,(*so.val_d)) ;
    double jlm2 = (ll<2) ? 0 : gsl_sf_bessel_jl (ll-2 ,(*so.val_d)) ;
    
    if (ll==0) {
      val = cos(*so.val_d)/(*so.val_d) - sin(*so.val_d)/(*so.val_d)/(*so.val_d) ;
    }
    else {
      val = jlm1 - (ll+1)/(*so.val_d)*jl ;
    }
    bool doder = (so.der_d == 0x0) ? false : true ;
    if (doder) {
      double der ;
      
      if (ll==0) { 
	der =  (-sin(*so.val_d)/(*so.val_d) -2*cos(*so.val_d)/(*so.val_d)/(*so.val_d) + 2*sin(*so.val_d)/(*so.val_d)/(*so.val_d)/(*so.val_d))*(*so.der_d) ;
      }
      else {
	if (ll==1) {
	der = (cos(*so.val_d)/(*so.val_d) -3*sin(*so.val_d)/(*so.val_d)/(*so.val_d) -6*cos(*so.val_d)/(*so.val_d)/(*so.val_d)/(*so.val_d) + 6*sin(*so.val_d) /(*so.val_d)/(*so.val_d)/(*so.val_d)/(*so.val_d))*(*so.der_d) ;
      }
	else {
	  der = (jlm2 - double(2*ll+1)*jlm1/(*so.val_d) + double((ll+1)*(ll+2))*jl/(*so.val_d)/(*so.val_d))*(*so.der_d) ;
	}
      }
      return Term_eq (dom, val, der) ;
    }
    else {
      return Term_eq (dom, val) ;
    }
  }
 
   // Scalar case
  if (so.type_data == TERM_T) {
    if (so.val_t->get_valence()!=0) {
      cerr << "Bessel function only defined for scalars so far..." << endl ;
      abort() ;
    }

    Val_domain val (so.val_t->get_space().get_domain(dom)) ;
    val.allocate_conf() ;
    Val_domain der (so.val_t->get_space().get_domain(dom)) ;
    der.allocate_conf() ;
  
    bool doder = (so.der_t == 0x0) ? false : true ;
    Index pos (val.get_domain()->get_nbr_points()) ;
    do  {
      double vv = (*so.val_t)()(dom)(pos) ;  
      
      double jl = gsl_sf_bessel_jl (ll, vv) ;
      double jlm1 = (ll==0) ? 0 : gsl_sf_bessel_jl (ll-1 , vv) ;
      double jlm2 = (ll<2) ? 0 : gsl_sf_bessel_jl (ll-2 , vv) ;
  
     if (ll==0) {
      val.set(pos) = cos(vv)/vv - sin(vv)/vv/vv ;
    }
    else {
      val.set(pos) = jlm1 - (ll+1)/vv*jl ;
    }
    
      if (doder) {
	  double dd = (*so.der_t)()(dom)(pos) ;
	   if (ll==0) { 
	der.set(pos) =  (-sin(vv)/vv -2*cos(vv)/vv/vv + 2*sin(vv)/vv/vv/vv)*dd ;
      }
      else {
	if (ll==1) {
	der.set(pos) = (cos(vv)/vv -3*sin(vv)/vv/vv -6*cos(vv)/vv/vv/vv + 6*sin(vv) /vv/vv/vv/vv)*dd ;
      }
	else {
	  der = (jlm2 - double(2*ll+1)*jlm1/vv + double((ll+1)*(ll+2))*jl/vv/vv)*dd ;
	}
	 
      }
      }
    }
    while (pos.inc()) ;
  
    val.std_base() ;
    if (doder)
      der.std_base() ;
  
    Scalar valres (so.val_t->get_space()) ;
    valres.set_domain(dom) = val ;
    if (doder) {
      Scalar derres (so.val_t->get_space()) ;
      derres.set_domain(dom) = der ;
      return Term_eq (dom, valres, derres) ;
    }
    else  {
      return Term_eq (dom, valres) ;
    }
  }
  
  // Unknown case
  cerr << "Unknown type of Term_eq in bessel_djl" << endl  ;
  abort() ;
}


Term_eq bessel_dyl (const Term_eq& so, int ll) {

  int dom = so.get_dom() ;
  
  // Double case 
  if (so.type_data == TERM_D) { 
    double val ;
    
    double yl = gsl_sf_bessel_yl (ll,(*so.val_d)) ;
    double ylm1 = (ll==0) ? 0 : gsl_sf_bessel_yl (ll-1 ,(*so.val_d)) ;
    double ylm2 = (ll<2) ? 0 : gsl_sf_bessel_yl (ll-2 ,(*so.val_d)) ;
    
    if (ll==0) {
      val = sin(*so.val_d)/(*so.val_d) + cos(*so.val_d)/(*so.val_d)/(*so.val_d) ;
    }
    else {
      val = ylm1 - (ll+1)/(*so.val_d)*yl ;
    }
    bool doder = (so.der_d == 0x0) ? false : true ;
    if (doder) {
      double der ;
      
      if (ll==0) { 
	der =  (cos(*so.val_d)/(*so.val_d) -2*sin(*so.val_d)/(*so.val_d)/(*so.val_d) - 2*cos(*so.val_d)/(*so.val_d)/(*so.val_d)/(*so.val_d))*(*so.der_d) ;
      }
      else {
	if (ll==1) {
	der = (sin(*so.val_d)/(*so.val_d) +3*cos(*so.val_d)/(*so.val_d)/(*so.val_d) -6*sin(*so.val_d)/(*so.val_d)/(*so.val_d)/(*so.val_d) -6*cos(*so.val_d) /(*so.val_d)/(*so.val_d)/(*so.val_d)/(*so.val_d))*(*so.der_d) ;
      }
	else {
	  der = (ylm2 - double(2*ll+1)*ylm1/(*so.val_d) + double((ll+1)*(ll+2))*yl/(*so.val_d)/(*so.val_d))*(*so.der_d) ;
	}
      }
      return Term_eq (dom, val, der) ;
    }
    else {
      return Term_eq (dom, val) ;
    }
  }
 
   // Scalar case
  if (so.type_data == TERM_T) {
    if (so.val_t->get_valence()!=0) {
      cerr << "Bessel function only defined for scalars so far..." << endl ;
      abort() ;
    }

    Val_domain val (so.val_t->get_space().get_domain(dom)) ;
    val.allocate_conf() ;
    Val_domain der (so.val_t->get_space().get_domain(dom)) ;
    der.allocate_conf() ;
  
    bool doder = (so.der_t == 0x0) ? false : true ;
    Index pos (val.get_domain()->get_nbr_points()) ;
    do  {
      double vv = (*so.val_t)()(dom)(pos) ;  
      
      double yl = gsl_sf_bessel_yl (ll, vv) ;
      double ylm1 = (ll==0) ? 0 : gsl_sf_bessel_yl (ll-1 , vv) ;
      double ylm2 = (ll<2) ? 0 : gsl_sf_bessel_yl (ll-2 , vv) ;
  
     if (ll==0) {
      val.set(pos) = sin(vv)/vv + cos(vv)/vv/vv ;
    }
    else {
      val.set(pos) = ylm1 - (ll+1)/vv*yl ;
    }
    
      if (doder) {
	  double dd = (*so.der_t)()(dom)(pos) ;
	   if (ll==0) { 
	der.set(pos) =  (cos(vv)/vv -2*sin(vv)/vv/vv - 2*cos(vv)/vv/vv/vv)*dd ;
      }
      else {
	if (ll==1) {
	der.set(pos) = (sin(vv)/vv +3*cos(vv)/vv/vv -6*sin(vv)/vv/vv/vv - 6*cos(vv) /vv/vv/vv/vv)*dd ;
      }
	else {
	  der = (ylm2 - double(2*ll+1)*ylm1/vv + double((ll+1)*(ll+2))*yl/vv/vv)*dd ;
	}
	 
      }
      }
    }
    while (pos.inc()) ;
  
    val.std_base() ;
    if (doder)
      der.std_base() ;
  
    Scalar valres (so.val_t->get_space()) ;
    valres.set_domain(dom) = val ;
    if (doder) {
      Scalar derres (so.val_t->get_space()) ;
      derres.set_domain(dom) = der ;
      return Term_eq (dom, valres, derres) ;
    }
    else  {
      return Term_eq (dom, valres) ;
    }
  }
  
  // Unknown case
  cerr << "Unknown type of Term_eq in bessel_dyl" << endl  ;
  abort() ;
}}
