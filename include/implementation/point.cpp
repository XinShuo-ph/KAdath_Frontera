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
#include "utilities.hpp"

namespace Kadath {
inline Point::Point(int dim) : ndim(dim) {
	coord = MemoryMapper::get_memory<double>(ndim);

	for (int i=0 ; i<ndim ; i++)
		coord[i] = 0 ;
}

inline Point::Point(const Point& so) : ndim(so.ndim) {
	coord = MemoryMapper::get_memory<double>(ndim);

	for (int i=0 ; i<ndim ; i++)
	   coord[i] = so.coord[i] ;
}

inline Point::Point (FILE* fd) {
	fread_be (&ndim, sizeof(int), 1, fd) ;

	coord = MemoryMapper::get_memory<double>(ndim);

	fread_be (coord, sizeof(double), ndim, fd) ;
}

inline Point::~Point() {
	MemoryMapper::release_memory<double>(coord,ndim);
}

inline void Point::save (FILE* fd) const {
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (coord, sizeof(double), ndim, fd) ;
}

inline void Point::operator= (const Point& so) {
	assert (ndim==so.ndim) ;
	for (int i=0 ; i<ndim ; i++)
	    coord[i] = so.coord[i] ;
}

inline double& Point::set(int i) {
	assert(i>0) ;
	assert(i<=ndim) ;
	return coord[i-1] ;
}

inline double Point::operator() (int i) const {
	assert(i>0) ;
	assert(i<=ndim) ;
	return coord[i-1] ;
}

}
