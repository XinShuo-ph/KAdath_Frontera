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

#include "adapted_polar.hpp"
#include "system_of_eqs.hpp"

namespace Kadath {

void Space_polar_adapted::add_eq (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused) const  {
	for (int dd=sys.get_dom_min() ; dd<sys.get_dom_max(); dd++) {
		sys.add_eq_inside (dd, eq, nused, pused) ;
		sys.add_eq_matching (dd, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (dd, OUTER_BC, rac_der, nused, pused) ;
	}
	sys.add_eq_inside (sys.get_dom_max(), eq, nused, pused) ;
}

}
