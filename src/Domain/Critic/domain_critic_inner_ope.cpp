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
#include "critic.hpp"
#include "array_math.hpp"
#include "val_domain.hpp"
namespace Kadath {
int div_x_1d (int, Array<double>&) ;

Val_domain Domain_critic_inner::der_partial_var (const Val_domain& so, int which_var) const {

	switch (which_var) {
		case 0 : 
			return so.der_var(1) / xlim ;
			break ;
		case 1 :
			return so.der_var(2) ;
			break ;
		
		default:
			cerr << "Unknown variable in Domain_critic_inner::der_partial_var" << endl ;
			abort() ;
		}
}

Val_domain Domain_critic_inner::div_x (const Val_domain& so) const {
	so.coef() ;
	Val_domain res(this) ;

	res.base= so.base ;
	
	res.cf = (so.base.ope_1d(div_x_1d, 0, so.cf, res.base)/xlim) ;

	res.in_coef = true ;
	return res ;
}}
