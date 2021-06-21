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
#include "polar.hpp"
#include "point.hpp"
#include "array.hpp"
#include "scalar.hpp"
#include "tensor.hpp"
namespace Kadath {
void Domain_polar_nucleus::affecte_tau_one_coef_val_domain (Val_domain& so, int mquant, int llim, int cc, int& conte) const {

	int lquant ;

	so.is_zero = false ;
	so.allocate_coef() ;
	*so.cf=0. ;
	Index pos_cf(nbr_coefs) ;

	bool found = false ;

		// Positions of the Galerkin basis
	Index pos_gal_t (nbr_coefs) ;
	Index pos_gal_r (nbr_coefs) ;
	Index pos_gal_rt (nbr_coefs) ;
	double fact_t, fact_r, fact_rt ;

	// Loop on theta
	int baset = (*so.get_base().bases_1d[1]) (0) ;
	for (int j=0 ; j<nbr_coefs(1) ; j++) {	
		int baser = (*so.get_base().bases_1d[0]) (j) ;
		pos_cf.set(1) = j ;
		// Loop on r :
		for (int i=0 ; i<nbr_coefs(0) ; i++) {
			pos_cf.set(0) = i ;
			switch (baset) {
				case COS_EVEN :
					lquant = 2*j ;
					// No galerkin :
					if ((mquant==0) && (lquant<=llim))  {
						if (conte==cc)  {
						  found = true ;
						  so.cf->set(pos_cf) = 1. ;
						  }
						  conte ++ ;
					}
					else if (mquant==0) {
					  if (i!=0) {
						if (conte==cc) {
						found = true ;
						// Galerkin base in r only
						pos_gal_r = pos_cf ;
						pos_gal_r.set(0) = 0 ;
						switch (baser) {
							case CHEB_EVEN :
							  fact_r = - pow(-1, i) ;
							  break ;
							case LEG_EVEN : {
							  fact_r = -1. ;
							  for (int t=0 ; t<i ; t++)
							    fact_r *= -double(2*t+1)/double(2*t+2) ;
							  }
							  break ;
							default :
							  cerr << "Strange base in Domain_polar_nucleus::affecte_one_coef_val_domain" << endl ;
							  abort()  ;
						}
						so.cf->set(pos_cf) = 1 ;
						so.cf->set(pos_gal_r) += fact_r ;
						}
							conte ++ ;
					  }
					}
					else if ((j!=0) && (i!=0)) {
					    if (conte==cc) {
					    found = true ;
					    // Need to use two_dimensional Galerkin basis (aouch !)
					    pos_gal_r = pos_cf ;
					    pos_gal_r.set(0) = 0 ;
					    pos_gal_t = pos_cf ;
					    pos_gal_t.set(1) = 0 ;
					    pos_gal_rt = pos_cf ;
					    pos_gal_rt.set(0) = 0 ;
					    pos_gal_rt.set(1) = 0 ;
					    switch (baser) {
						case CHEB_EVEN :
							fact_r = -pow(-1, i) ;
							fact_t = -1. ;
							fact_rt = pow(-1, i) ;
							break ;
						  case LEG_EVEN : {
							double l0 = 1 ;
							 for (int t=0 ; t<i ; t++)
							    l0 *= -double(2*t+1)/double(2*t+2) ;
							fact_r = - l0 ;
							fact_t = -1. ;
							fact_rt = l0 ;
							}
						break ;
						 default :
							  cerr << "Strange base in Domain_polar_nucleus::affecte_one_coef_val_domain" << endl ;
							  abort()  ;
						}
						  so.cf->set(pos_cf) = 1. ;
						  so.cf->set(pos_gal_r) = fact_r ;
						  so.cf->set(pos_gal_t) = fact_t ;
						  so.cf->set(pos_gal_rt) = fact_rt ;
					      }
					   conte ++ ;
					   }
					break ;
			case COS_ODD:
				lquant = 2*j+1 ;
				if ((j!=nbr_coefs(1)-1) && (i!=nbr_coefs(0)-1))  {
					      if ((mquant==0) && (lquant<=llim+1)) {
						  if (conte==cc) {
						      found = true ;
						      so.cf->set(pos_cf) = 1. ;
						     }
						  conte++  ;
						  }
						 else {
							if ((mquant==0) && (i!=0)) {
							 if (conte==cc) {
							  found = true ;
							   pos_gal_r = pos_cf ;
							  pos_gal_r.set(0) = 0 ;
							  switch (baser) {
							    case CHEB_ODD :
							      fact_r = - (2*i+1) * pow(-1, i) ;
								break ;
							  case LEG_ODD : {
								fact_r = -1. ;
								for (int t=0 ; t<i ; t++)
								  fact_r *= -double(2*t+3)/double(2*t+2) ;
								}
								break ;
							  default :
								cerr << "Strange base in Domain_polar_nucleus::affecte_one_coef_val_domain" << endl ;
								abort()  ;
							   }

						  so.cf->set(pos_cf) = 1. ;
						  so.cf->set(pos_gal_r) = fact_r ;
						}
						   conte ++ ;
					}
					else if ((j!=0) && (i!=0)) {
					  if (conte==cc) {
						found = true ;
						 // Need to use two_dimensional Galerkin basis (aouch !)
						pos_gal_r = pos_cf ;
						pos_gal_r.set(0) = 0 ;
						pos_gal_t = pos_cf ;
						pos_gal_t.set(1) = 0 ;
						pos_gal_rt = pos_cf ;
						pos_gal_rt.set(0) = 0 ;
						pos_gal_rt.set(1) = 0 ;
						switch (baser) {
						  case CHEB_ODD :
						    fact_r = -pow(-1, i)*(2*i+1) ;
						    fact_t = -1. ;
						    fact_rt = pow(-1, i)*(2*i+1) ;
						    break ;
						  case LEG_ODD : {
						    double l0 = 1 ;
						    for (int t=0 ; t<i ; t++)
						      l0 *= -double(2*t+3)/double(2*t+2) ;
						    fact_r = - l0 ;
						    fact_t = -1. ;
						    fact_rt = l0 ;
						    }
						    break ;
						default :
							  cerr << "Strange base in Domain_polar_nucleus::affecte_one_coef_val_domain" << endl ;
							  abort()  ;
						}  
						 so.cf->set(pos_cf) = 1. ;
						 so.cf->set(pos_gal_r) = fact_r ;
					         so.cf->set(pos_gal_t) = fact_t ;
						 so.cf->set(pos_gal_rt) = fact_rt ;
						  }
						  conte ++ ;
						}
					}
				}
				break ;
			case SIN_EVEN:
				lquant = 2*j ;
				if ((j!=0) && (j!=nbr_coefs(1)-1)) { 
						if ((mquant<=1) && (lquant<=llim)) {
						    if (conte==cc) {
								      found = true ;
								      so.cf->set(pos_cf) = 1. ;
								      }
						    conte ++ ;
						}
						else {
						if ((mquant<=1) && (i!=0)) {
							// Galerkin base in r only
							 if (conte==cc) {
								found = true ;
							pos_gal_r = pos_cf ;
							pos_gal_r.set(0) = 0 ;
							switch (baser) {
								case CHEB_EVEN :
								  fact_r = - pow(-1, i) ;
								  break ;
								case LEG_EVEN : {
								  fact_r = -1. ;
								  for (int t=0 ; t<i ; t++)
								  fact_r *= -double(2*t+1)/double(2*t+2) ;
								  }
								  break ;
								default :
								  cerr << "Strange base in Domain_polar_nucleus::affecte_tau_val_domain" << endl ;
								  abort()  ;
								}	  
							    so.cf->set(pos_cf) = 1. ;
								so.cf->set(pos_gal_r) = fact_r ;
								}
							    conte ++ ;
							}
						else  {
						//Double Galerkin
						if ((j!=1) && (i!=0)) {

								 if (conte==cc) {
								    found = true ;
								 // Need to use two_dimensional Galerkin basis (aouch !)
							    pos_gal_r = pos_cf ;
							    pos_gal_r.set(0) = 0 ;
							    pos_gal_t = pos_cf ;
							    pos_gal_t.set(1) = 1 ;
							    pos_gal_rt = pos_cf ;
							    pos_gal_rt.set(0) = 0 ;
							    pos_gal_rt.set(1) = 1 ;
							     switch (baser) {
							      case CHEB_EVEN :
								fact_r = -pow(-1, i) ;
								fact_t = -j ;
								fact_rt = pow(-1, i)*j ;
								break ;
							      case LEG_EVEN : {
								double l0 = 1 ;
								 for (int t=0 ; t<i ; t++)
								  l0 *= -double(2*t+1)/double(2*t+2) ;
								fact_r = - l0 ;
								fact_t = -j ;
								fact_rt = l0*j ;
								}
								break ;
							    default :
									  cerr << "Strange base in Domain_polar_nucleus::affecte_tau_val_domain" << endl ;
									  abort()  ;
								}  
							      so.cf->set(pos_cf) = 1 ;
							      so.cf->set(pos_gal_r) = fact_r ;
							      so.cf->set(pos_gal_t) = fact_t ;
							      so.cf->set(pos_gal_rt) = fact_rt ;
								}
							      conte ++ ;
							}
						      }
						}
						}
				break ;
			case SIN_ODD:
				lquant = 2*j+1 ;
				if ((j!=nbr_coefs(1)-1) && (i!=nbr_coefs(0)-1))  {
							      if ((mquant<=1) && (lquant<=llim+1)) {
								  if (conte==cc) {
								      found = true ;
								      so.cf->set(pos_cf) = 1. ;
								      }
								  conte++  ;
							      }
							      else {
								if ((mquant<=1) && (i!=0)) {
								  if (conte==cc) {
								    found = true ;
								    pos_gal_r = pos_cf ;
								    pos_gal_r.set(0) = 0 ;
								    switch (baser) {
								      case CHEB_ODD :
									fact_r = - (2*i+1) * pow(-1, i) ;
									break ;
								      case LEG_ODD : {
									fact_r = -1. ;
									for (int t=0 ; t<i ; t++)
									  fact_r *= -double(2*t+3)/double(2*t+2) ;
									}
									break ;
								      default :
									cerr << "Strange base in Domain_polar_nucleus::affecte_one_coef_val_domain" << endl ;
									abort()  ;
								      }

								  so.cf->set(pos_cf) = 1. ;
								  so.cf->set(pos_gal_r) = fact_r ;
								}
							       conte ++ ;
							}
							else if ((j!=0) && (i!=0)) {
							    if (conte==cc) {
								found = true ;
								 // Need to use two_dimensional Galerkin basis (aouch !)
								pos_gal_r = pos_cf ;
								pos_gal_r.set(0) = 0 ;
								pos_gal_t = pos_cf ;
								pos_gal_t.set(1) = 0 ;
								pos_gal_rt = pos_cf ;
								pos_gal_rt.set(0) = 0 ;
								pos_gal_rt.set(1) = 0 ;
								switch (baser) {
								  case CHEB_ODD :
								    fact_r = -pow(-1, i)*(2*i+1) ;
								    fact_t = -(2*j+1) ;
								    fact_rt = pow(-1, i)*(2*i+1)*(2*j+1) ;
								    break ;
								  case LEG_ODD : {
								    double l0 = 1 ;
								    for (int t=0 ; t<i ; t++)
								      l0 *= -double(2*t+3)/double(2*t+2) ;
								    fact_r = - l0 ;
								    fact_t = -(2*j+1) ;
								    fact_rt = l0*(2*j+1) ;
								    }
								    break ;
								default :
									  cerr << "Strange base in Domain_polar_nucleus::affecte_one_coef_val_domain" << endl ;
									  abort()  ;
								}  
							      so.cf->set(pos_cf) = 1. ;
							      so.cf->set(pos_gal_r) = fact_r ;
							      so.cf->set(pos_gal_t) = fact_t ;
							      so.cf->set(pos_gal_rt) = fact_rt ;
							      }
							      conte ++ ;
							}
						    }
						}
				break ;
			default:
				cerr << "Unknow theta basis in Domain_polar_nucleus::affecte_coef_val_domain" << endl ;
				abort() ;
			}
			}
	}
	// If not found put to zero :
	if (!found)
		so.set_zero() ;
}

void Domain_polar_nucleus::affecte_tau_one_coef (Tensor& tt, int dom, int cc, int& pos_cf) const {

	// Check right domain
	assert (tt.get_space().get_domain(dom)==this) ;

	int val = tt.get_valence() ;
	switch (val) {
		case 0 :
			if (!tt.is_m_quant_affected())
			  affecte_tau_one_coef_val_domain (tt.set().set_domain(dom), 0, 0, cc, pos_cf) ;
			else
			  affecte_tau_one_coef_val_domain (tt.set().set_domain(dom), tt.get_parameters()->get_m_quant(), 0, cc, pos_cf) ;
			break ;
	
		default :
			cerr << "Valence " << val << " not implemented in Domain_polar_nucleus::affecte_tau" << endl ;
			break ;
	}
}}
