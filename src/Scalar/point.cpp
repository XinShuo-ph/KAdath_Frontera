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
#include <assert.h>
#include "point.hpp"
#include "utilities.hpp"
namespace Kadath {

ostream& operator<< (ostream& o, const Point& so) {

	o << "(" ;
	for (int i=0 ; i<so.ndim-1 ; i++)
	    o << so.coord[i] << ", " ;
	o << so.coord[so.ndim-1] << ")" ;
        return o ;
}}
