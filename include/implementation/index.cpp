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

#include "assert.h"
#include "dim_array.hpp"

namespace Kadath {
// Constructor
inline Index::Index(const Dim_array& res) : sizes(res) {
  coord = MemoryMapper::get_memory<int>(get_ndim());

	for (int i=0 ; i<get_ndim() ; i++)
	    coord[i] = 0 ;
}

// Copy constructor
inline Index::Index(const Index& so) : sizes(so.sizes) {
  coord = MemoryMapper::get_memory<int>(get_ndim());

	for (int i=0 ; i<get_ndim() ; i++)
	   coord[i] = so.coord[i] ;
}

// Destructor
inline Index::~Index() {
	MemoryMapper::release_memory<int>(coord,get_ndim());
}

// Read/write
inline int& Index::set(int i) {
	assert(i>=0) ;
	assert(i<get_ndim()) ;
	return coord[i] ;
}

// Read only
inline int Index::operator() (int i) const {
	assert(i>=0) ;
	assert(i<get_ndim()) ;
	return coord[i] ;
}

// Assignement
inline void Index::operator= (const Index& so) {
	assert (sizes == so.sizes) ;
	for (int i=0 ; i<get_ndim() ; i++)
	    coord[i] = so.coord[i] ;
}

inline void Index::set_start() {
	for (int i=0 ; i<get_ndim() ; i++)
	     coord[i] = 0 ;
}

// Increment
inline bool Index::inc (int increm, int var) {
	bool res = ((var >=0) && (var<get_ndim())) ? true : false ;
	if (res) {
	coord[var] += increm ;
	for (int i=var ; i<get_ndim()-1 ; i++) {
		div_t division = div(coord[i], sizes(i)) ;
		coord[i] = division.rem ;
		coord[i+1] += division.quot ;
	}
        res = (coord[get_ndim()-1] >= sizes(get_ndim()-1)) ? false : true ;
	}
	return res ;
}

inline bool Index::inc () {
  for(int i = 0; i < get_ndim(); ++i) {
    if (coord[i]+1 < sizes(i)) {
      coord[i]++;
      return true;
    }
    else {
      coord[i] = 0;
    }
  }
  return false;
}


inline bool Index::operator== (const Index& xx) const {
	bool res = (get_ndim()==xx.get_ndim()) ? true : false ;
	if (res)
		for (int i=0 ; i<get_ndim() ; i++)
			if (xx.coord[i] != coord[i])
				res = false ;
	return res ;
}

}
