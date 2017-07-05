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

namespace Kadath {
Term_eq operator+ (const Term_eq& aa, double xx) {
	Term_eq auxi (aa.dom, xx, 0.) ;
	return aa + auxi ;
}

Term_eq operator+ (double xx, const Term_eq& aa) {
	Term_eq auxi (aa.dom, xx, 0.) ;
	return aa + auxi ;
}

Term_eq operator- (const Term_eq& aa) {
	return (-1. * aa) ;
}

Term_eq Term_eq::der_abs (int i) const {

	Term_eq res (dom, type_data) ;
	Index pos(*val_t) ;
	switch (type_data) {
		case TERM_D :
			cerr << "derivative only defined with respect to tensors" << endl ;
			abort() ;
		case TERM_T :
			res.val_t = new Tensor(*val_t, false) ;
			if (der_t!=0x0)
				res.der_t = new Tensor(*der_t, false) ;
			do {
				res.val_t->set(pos).set_domain(dom) = (*val_t)(pos)(dom).der_abs(i) ;
				if (der_t !=0x0)
					res.der_t->set(pos).set_domain(dom) = (*der_t)(pos)(dom).der_abs(i) ;
			}
			while (pos.inc()) ;
			break ;
		default: 
			cerr << "Unknown data storage in operator+" << endl ;
			abort() ;
		}	
	return res ;
}

Term_eq operator+ (const Term_eq& aa, const Term_eq& bb) {
	

	bool do_der = true ;

	int data_res ;
	assert (aa.dom==bb.dom) ;
	switch (aa.type_data) {
		case (TERM_D) :
			if (aa.der_d==0x0)
				do_der = false ;
			switch (bb.type_data) {
				case (TERM_D) :
					data_res = TERM_D ;
					if (bb.der_d==0x0)
						do_der=false ;
					break ;
				case (TERM_T) :
					data_res = TERM_T ;
					if (bb.der_t==0x0)
						do_der = false ;
					break ;
				default: 
					cerr << "Unknown data storage in operator+" << endl ;
					abort() ;
				}
			break ;
		case (TERM_T) :
			if (aa.der_t==0x0)
				do_der = false ;
			switch (bb.type_data) {
				case (TERM_D) :
					if (bb.der_d==0x0)
						do_der=false ;
					data_res = TERM_T ;
					break ;
				case (TERM_T) :
					if (bb.der_t==0x0)
						do_der=false ;
					data_res = TERM_T ;
					break ;
				default: 
					cerr << "Unknown data storage in operator+" << endl ;
					abort() ;
				}
			break ;
		default :
			cerr << "Unknown data storage in operator+" << endl ;
			abort() ;	
	}

	Term_eq res(aa.dom, data_res) ;
	switch (aa.type_data) {
		case (TERM_D) :
			switch (bb.type_data) {
				case (TERM_D) :
					res.val_d = new double(*aa.val_d + *bb.val_d) ;
					res.der_d = (do_der) ? new double (*aa.der_d + *bb.der_d) :0x0 ;
					break ;
				case (TERM_T) :
					res.val_t = new Tensor(add_one_dom(res.dom, *aa.val_d, *bb.val_t)) ;
					res.der_t = (do_der) ? 
						new Tensor (add_one_dom(res.dom, *aa.der_d, *bb.der_t)) :0x0;
					break ;
				default: 
					cerr << "Unknown data storage in operator+" << endl ;
					abort() ;
				}
			break ;
		case (TERM_T) :
			switch (bb.type_data) {
				case (TERM_D) :
					res.val_t = new Tensor(add_one_dom(res.dom, *aa.val_t , *bb.val_d)) ;
					res.der_t = (do_der) ? 
						new Tensor(add_one_dom(res.dom, *aa.der_t , *bb.der_d)) :0x0;
					break ;
				case (TERM_T) :
					res.val_t = new Tensor(add_one_dom(res.dom, *aa.val_t , *bb.val_t)) ;
					res.der_t = (do_der) ? 
					 new Tensor(add_one_dom(res.dom, *aa.der_t , *bb.der_t)) : 0x0 ;
					break ;
				default: 
					cerr << "Unknown data storage in operator+" << endl ;
					abort() ;
				}
			break ;
		default :
			cerr << "Unknown data storage in operator+" << endl ;
			abort() ;	
	}
	return res ;
}


Term_eq operator- (const Term_eq& aa, const Term_eq& bb) {
	assert (aa.dom==bb.dom) ;
	int data_res ;
	bool do_der = true ;
	switch (aa.type_data) {
		case (TERM_D) :
			if (aa.der_d==0x0)
				do_der = false ;
			switch (bb.type_data) {
				case (TERM_D) :
					if (bb.der_d==0x0)
						do_der = false ;
					data_res = TERM_D ;
					break ;
				case (TERM_T) :
					if (bb.der_t==0x0)
						do_der = false ;
					data_res = TERM_T ;
					break ;
				default: 
					cerr << "Unknown data storage in operator-" << endl ;
					abort() ;
				}
			break ;
		case (TERM_T) :
			if (aa.der_t==0x0)
				do_der = false ;
			switch (bb.type_data) {
				case (TERM_D) :
					if (bb.der_d==0x0)
						do_der = false ;
					data_res = TERM_T ;
					break ;
				case (TERM_T) :
					if (bb.der_t==0x0)
						do_der = false ;
					data_res = TERM_T ;
					break ;
				default: 
					cerr << "Unknown data storage in operator-" << endl ;
					abort() ;
				}
			break ;
		default :
			cerr << "Unknown data storage in operator-" << endl ;
			abort() ;	
	}

	Term_eq res(aa.dom, data_res) ;
	switch (aa.type_data) {
		case (TERM_D) :
			switch (bb.type_data) {
				case (TERM_D) :
					res.val_d = new double(*aa.val_d - *bb.val_d) ;
					res.der_d = (do_der) ? new double (*aa.der_d - *bb.der_d) :0x0 ;
					break ;
				case (TERM_T) :
					res.val_t = new Tensor(sub_one_dom(res.dom, *aa.val_d, *bb.val_t)) ;
					res.der_t = (do_der) ? 
						new Tensor (sub_one_dom(res.dom, *aa.der_d, *bb.der_t)) : 0x0 ;
					break ;
				default: 
					cerr << "Unknown data storage in operator+" << endl ;
					abort() ;
				}
			break ;
		case (TERM_T) :
			switch (bb.type_data) {
				case (TERM_D) :
					res.val_t = new Tensor(sub_one_dom(res.dom, *aa.val_t, *bb.val_d)) ;
					res.der_t = (do_der) ? 
						new Tensor(sub_one_dom(res.dom, *aa.der_t, *bb.der_d)) : 0x0 ;
					break ;
				case (TERM_T) :
					res.val_t = new Tensor(sub_one_dom(res.dom, *aa.val_t, *bb.val_t)) ;
					res.der_t = (do_der) ? 
						new Tensor(sub_one_dom(res.dom, *aa.der_t, *bb.der_t)) : 0x0 ;
					break ;
				default: 
					cerr << "Unknown data storage in operator+" << endl ;
					abort() ;
				}
			break ;
		default :
			cerr << "Unknown data storage in operator+" << endl ;
			abort() ;	
	}
	return res ;
}

Term_eq operator* (const Term_eq& aa, const Term_eq& bb) {
	
	assert (aa.dom==bb.dom) ;
	int data_res ;
	bool do_der = true ;
	switch (aa.type_data) {
		case (TERM_D) :
			if (aa.der_d==0x0)
				do_der = false ;
			switch (bb.type_data) {
				case (TERM_D) :
					if (bb.der_d==0x0)
						do_der = false ;
					data_res = TERM_D ;
					break ;
				case (TERM_T) :
					if (bb.der_t==0x0)
						do_der = false ;
					data_res = TERM_T ;
					break ;
				default: 
					cerr << "Unknown data storage in operator*" << endl ;
					abort() ;
				}
			break ;
		case (TERM_T) :
			if (aa.der_t==0x0)
				do_der=false ;
			switch (bb.type_data) {
				case (TERM_D) :
					if (bb.der_d==0x0)
						do_der=false ;
					data_res = TERM_T ;
					break ;
				case (TERM_T) :
					if (bb.der_t==0x0)
						do_der = false ;
					data_res = TERM_T ;
					break ;
				default: 
					cerr << "Unknown data storage in operator*" << endl ;
					abort() ;
				}
			break ;
		default :
			cerr << "Unknown data storage in operator*" << endl ;
			abort() ;	
	}

	Term_eq res(aa.dom, data_res) ;
		switch (aa.type_data) {
		case (TERM_D) :
			switch (bb.type_data) {
				case (TERM_D) :
					res.val_d = new double((*aa.val_d) * (*bb.val_d)) ;
					res.der_d = (do_der) ? 
						new double ((*aa.der_d)*(*bb.val_d) + (*aa.val_d)*(*bb.der_d)) : 0x0 ;
					break ;
				case (TERM_T) :
					res.val_t = new Tensor(mult_one_dom(res.dom, *aa.val_d, *bb.val_t)) ;
					res.der_t = (do_der) ? 
	new Tensor (add_one_dom(res.dom, mult_one_dom(res.dom, *aa.der_d,*bb.val_t), 
				mult_one_dom(res.dom, *aa.val_d, *bb.der_t))) : 0x0 ;
					break ;
				default: 
					cerr << "Unknown data storage in operator+" << endl ;
					abort() ;
				}
			break ;
		case (TERM_T) :
			switch (bb.type_data) {
				case (TERM_D) :
					res.val_t = new Tensor(mult_one_dom(res.dom, *aa.val_t, *bb.val_d)) ;
					res.der_t = (do_der) ? 
new Tensor(add_one_dom(res.dom, mult_one_dom(res.dom, *aa.der_t, *bb.val_d), mult_one_dom(res.dom, *aa.val_t, *bb.der_d))) : 0x0 ;
					break ;
				case (TERM_T) :
					res.val_t = new Tensor(mult_one_dom(res.dom, *aa.val_t, *bb.val_t)) ;
					res.der_t = (do_der) ? 
new Tensor(add_one_dom(res.dom, 
	mult_one_dom(res.dom, *aa.der_t, *bb.val_t), mult_one_dom(res.dom, *aa.val_t, *bb.der_t))) : 0x0 ;
					break ;
				default: 
					cerr << "Unknown data storage in operator*" << endl ;
					abort() ;
				}
			break ;
		default :
			cerr << "Unknown data storage in operator*" << endl ;
			abort() ;	
	}

	return res ;
}


Term_eq operator/ (const Term_eq& aa, const Term_eq& bb) {
	
	assert (aa.dom==bb.dom) ;
	int data_res ;
	bool do_der = true ;
	switch (aa.type_data) {
		case (TERM_D) :
			if (aa.der_d==0x0)
				do_der = false ;
			switch (bb.type_data) {
				case (TERM_D) :
					if (bb.der_d==0x0)
						do_der = false ;
					data_res = TERM_D ;
					break ;
				case (TERM_T) :
					if (bb.der_t==0x0)
						do_der = false ;
					data_res = TERM_T ;
					break ;
				default: 
					cerr << "Unknown data storage in operator*" << endl ;
					abort() ;
				}
			break ;
		case (TERM_T) :
			if (aa.der_t==0x0)
				do_der=false ;
			switch (bb.type_data) {
				case (TERM_D) :
					if (bb.der_d==0x0)
						do_der = false ;
					data_res = TERM_T ;
					break ;
				case (TERM_T) :
					if (bb.der_t==0x0)
						do_der = false ;
					data_res = TERM_T ;
					break ;
				default: 
					cerr << "Unknown data storage in operator*" << endl ;
					abort() ;
				}
			break ;
		default :
			cerr << "Unknown data storage in operator*" << endl ;
			abort() ;	
	}

	Term_eq res(aa.dom, data_res) ;
	switch (aa.type_data) {
		case (TERM_D) :
			switch (bb.type_data) {
				case (TERM_D) :
					res.val_d = new double((*aa.val_d) / (*bb.val_d)) ;
					res.der_d = (do_der) ? 
	new double (((*aa.der_d)*(*bb.val_d) - (*aa.val_d)*(*bb.der_d))/
								((*bb.val_d)*(*bb.val_d))) : 0x0 ;
					break ;
				case (TERM_T) :
					res.val_t = new Tensor(div_one_dom(res.dom, *aa.val_d, *bb.val_t)) ;
					res.der_t = (do_der) ? new Tensor (div_one_dom(res.dom, 
	sub_one_dom(res.dom, mult_one_dom(res.dom, *aa.der_d, *bb.val_t), 
			     mult_one_dom(res.dom, *aa.val_d, *bb.der_t)), 
	mult_one_dom(res.dom, *bb.val_t, *bb.val_t))) : 0x0 ;
					break ;
				default: 
					cerr << "Unknown data storage in operator+" << endl ;
					abort() ;
				}
			break ;
		case (TERM_T) :
			switch (bb.type_data) {
				case (TERM_D) :
					res.val_t = new Tensor(div_one_dom(res.dom, *aa.val_t, *bb.val_d)) ;
					res.der_t = (do_der) ? new Tensor(div_one_dom(res.dom, sub_one_dom(res.dom, 
			     mult_one_dom(res.dom, *aa.der_t, *bb.val_d), 
			     mult_one_dom(res.dom, *aa.val_t, *bb.der_d)), (*bb.val_d) * (*bb.val_d))) :  0x0 ;
					break ;
				case (TERM_T) :
					res.val_t = new Tensor(div_one_dom(res.dom, *aa.val_t, *bb.val_t)) ;
					res.der_t = (do_der) ? new Tensor(div_one_dom(res.dom, 
	sub_one_dom(res.dom, mult_one_dom(res.dom, *aa.der_t, *bb.val_t), 
			     mult_one_dom(res.dom, *aa.val_t, *bb.der_t)), 
	mult_one_dom(res.dom, *bb.val_t, *bb.val_t))) : 0x0 ;
					break ;
				default: 
					cerr << "Unknown data storage in operator+" << endl ;
					abort() ;
				}
			break ;
		default :
			cerr << "Unknown data storage in operator+" << endl ;
			abort() ;	
	}
	return res ;
}


Term_eq scalar_product (const Term_eq& aa, const Term_eq& bb) {
	
	assert (aa.dom==bb.dom) ;
	bool do_der = true ;
	switch (aa.type_data) {
		case (TERM_D) :
			cerr << "scalar_product only defined with tensors" << endl ;
			abort() ;
			break ;
		case (TERM_T) :
			if (aa.der_t==0x0)
				do_der=false ;
			switch (bb.type_data) {
				case (TERM_D) :
					cerr << "scalar_product only defined with tensors" << endl ;
					abort() ;
					break ;
				case (TERM_T) :
					if (bb.der_t==0x0)
						do_der = false ;
					break ;
				default: 
					cerr << "Unknown data storage in operator*" << endl ;
					abort() ;
				}
			break ;
		default :
			cerr << "Unknown data storage in operator*" << endl ;
			abort() ;	
	}

	// Try to see the parameter :
	
	bool param_val = false ;
	int mval = 0 ;
	
	bool param_der = false ;
	int mder = 0 ;
	if (aa.type_data==TERM_T) {
	  
	  if (aa.val_t->is_m_quant_affected())
	    param_val = true ;
	  if (bb.val_t->is_m_quant_affected())
	    param_val = true ;
	  if (param_val) {
	      if (aa.val_t->is_m_quant_affected())
		mval += aa.val_t->get_parameters()->get_m_quant() ;
	      if (bb.val_t->is_m_quant_affected())
		mval += bb.val_t->get_parameters()->get_m_quant() ;
	  }
	  
	  if (do_der) {
	      if (aa.der_t->is_m_quant_affected())
		param_der = true ;
	      if (bb.der_t->is_m_quant_affected())
		param_der = true ;
	  if (param_der) {
	      if (aa.der_t->is_m_quant_affected())
		mder += aa.der_t->get_parameters()->get_m_quant() ;
	      if (bb.der_t->is_m_quant_affected())
		mder += bb.der_t->get_parameters()->get_m_quant() ;
	  }
	  }
	}

	Term_eq res(aa.dom, TERM_T) ;
	switch (aa.type_data) {
		case (TERM_D) :
			cerr << "scalar_product only defined with tensors" << endl ;
			abort() ;
			break ;
		case (TERM_T) :
			switch (bb.type_data) {
				case (TERM_D) :
					cerr << "scalar_product only defined with tensors" << endl ;
					abort() ;
					break ;
				case (TERM_T) :
					res.val_t = new Tensor(scal_one_dom(res.dom, *aa.val_t, *bb.val_t)) ;
					if (mval != 0) {
					    res.val_t->affect_parameters() ;
					    res.val_t->set_parameters()->set_m_quant() = mval ;
					}
					
					res.der_t = (do_der) ? new Tensor(add_one_dom(res.dom, 
		scal_one_dom(res.dom, *aa.der_t, *bb.val_t), scal_one_dom(res.dom, *aa.val_t, *bb.der_t))) : 0x0 ;
		
					if (mder != 0) {
					    res.der_t->affect_parameters() ;
					    res.der_t->set_parameters()->set_m_quant() = mder ;
					}
					break ;
				default: 
					cerr << "Unknown data storage in operator+" << endl ;
					abort() ;
				}
			break ;
		default :
			cerr << "Unknown data storage in operator+" << endl ;
			abort() ;	
	}
	return res ;
}

Term_eq pow (const Term_eq& so, int nn) {
	
	if (nn>0) {
		Term_eq res (so) ;
		for (int i=1 ; i<nn ; i++)
			res = res*so ;
		return res ;
	}
	else {
		Scalar one_scal (so.get_val_t(), false) ;
		one_scal = 1. ;
		one_scal.std_base() ;
		Scalar zero_scal (one_scal) ;
		zero_scal = 0. ;
		Term_eq one (so.get_dom(), one_scal, zero_scal) ;
		Term_eq res (one) ;
		for (int i=0 ; i<-nn ; i++)
			res = res / so ;
		return res ;
	}
}

Term_eq operator* (int nn, const Term_eq& so) {
	
	Term_eq res (so.dom, so.type_data) ;

	switch (so.type_data) {
		case (TERM_D) :
			res.val_d = new double(double(nn)*(*so.val_d)) ;
			if (so.der_d!=0x0)
				res.der_d = new double(double(nn)*(*so.der_d)) ;
			break ;
		case (TERM_T) :
			res.val_t = new Tensor(mult_one_dom(so.dom, nn, *so.val_t)) ;
			if (so.der_t!=0x0)
				res.der_t = new Tensor(mult_one_dom(so.dom, nn, (*so.der_t))) ;
			break ;
		default :
			cerr << "Unknown data storage in operator+" << endl ;
			abort() ;
		}
	return res ;
}

Term_eq operator* (const Term_eq& so, int nn) {
	return nn*so ;
}

Term_eq operator* (double xx, const Term_eq& so) {
	
	Term_eq res (so.dom, so.type_data) ;

	switch (so.type_data) {
		case (TERM_D) :
			res.val_d = new double(xx*(*so.val_d)) ;
			if (so.der_d!=0x0)
				res.der_d = new double(xx*(*so.der_d)) ;
			break ;
		case (TERM_T) :
			res.val_t = new Tensor(mult_one_dom(so.dom, xx, (*so.val_t))) ;
			if (so.der_t!=0x0)
				res.der_t = new Tensor(mult_one_dom(so.dom, xx, (*so.der_t))) ;
			break ;
		default :
			cerr << "Unknown data storage in operator+" << endl ;
			abort() ;
		}
	return res ;
}

Term_eq operator* (const Term_eq& so, double xx) {
	return xx*so ;
}

Term_eq operator/ (const Term_eq& so, double xx) {
	
	Term_eq res (so.dom, so.type_data) ;

	switch (so.type_data) {
		case (TERM_D) :
			res.val_d = new double((*so.val_d)/xx) ;
			if (so.der_d!=0x0)
				res.der_d = new double((*so.der_d)/xx) ;
			break ;
		case (TERM_T) :
			res.val_t = new Tensor(div_one_dom(so.dom, (*so.val_t), xx)) ;
			if (so.der_t!=0x0)
				res.der_t = new Tensor(div_one_dom(so.dom, (*so.der_t), xx)) ;
			break ;
		default :
			cerr << "Unknown data storage in operator/" << endl ;
			abort() ;
		}
	return res ;
}

Term_eq sqrt (const Term_eq& so) {
	
	Term_eq res (so.dom, so.type_data) ;

	switch (so.type_data) {
		case (TERM_D) :
			res.val_d = new double(sqrt(*so.val_d)) ;
			if (so.der_d!=0x0)
				res.der_d = new double((*so.der_d)/2./sqrt(*so.val_d)) ;
			break ;
		case (TERM_T) :
			res.val_t = new Tensor(sqrt_one_dom(so.dom, *so.val_t)) ;
			if (so.der_t!=0x0)
				res.der_t = new Tensor(div_one_dom(so.dom, (*so.der_t), 2*sqrt_one_dom(so.dom, *so.val_t))) ;
			break ;
		default :
			cerr << "Unknown data storage in operator/" << endl ;
			abort() ;
		}
	return res ;
}

Term_eq partial (const Term_eq& so, char ind) {
	Term_eq res (so.dom, so.type_data) ;
	switch (so.type_data) {
		case (TERM_D) :
			cerr << "Partial only defined with respect to tensors" << endl ;
			abort() ;
		case (TERM_T) :
			res.val_t = new Tensor(partial_one_dom(so.dom, ind, *so.val_t)) ;
			if (so.der_t!=0x0)
				res.der_t = new Tensor(partial_one_dom(so.dom, ind, (*so.der_t))) ;
			break ;
		default :
			cerr << "Unknown data storage in partial" << endl ;
			abort() ;
		}
	return res ;
}
}
