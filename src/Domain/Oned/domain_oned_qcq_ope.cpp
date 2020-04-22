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
#include "oned.hpp"
#include "array_math.hpp"
#include "val_domain.hpp"
namespace Kadath {

int div_xp1_1d (int, Array<double>&) ;

Val_domain Domain_oned_qcq::der_partial_var (const Val_domain& so, int which_var) const {

	switch (which_var) {
		case 0 : 
			return so.der_var(1) / alpha ;
			break ;
		default:
			cerr << "Unknown variable in Domain_oned_qcq::der_partial_var" << endl ;
			abort() ;
		}
}

Val_domain Domain_oned_qcq::div_xp1 (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;

	res.cf = new Array<double> (so.base.ope_1d(div_xp1_1d, 0, *so.cf, res.base)) ;
	res.in_coef = true ;
	return res ;
}

double Domain_oned_qcq::integrale (const Val_domain& so) const {
    double res = 0 ;
    so.coef() ;
    int baser = (*so.get_base().bases_1d[0]) (0) ;
    
    switch (baser) {
      case CHEB :
	for (int i=0 ; i<nbr_coefs(0) ; i+=2)
	    res += so.get_coef()(i)*(1./double(i+1) - 1./double(i-1)) ;
	break ;
      default : 
	 cerr << "Case not implemented in Domain_oned_qcq::integrale" << endl ;
    }
    return res*alpha ;
}}
