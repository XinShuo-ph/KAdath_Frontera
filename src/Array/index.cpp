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
#include "index.hpp"
#include "dim_array.hpp"
#include "tensor.hpp"
namespace Kadath {
// Constructor
Index::Index(const Dim_array& res) : sizes(res) {
	coord = new int[get_ndim()] ;
	for (int i=0 ; i<get_ndim() ; i++)
	    coord[i] = 0 ;
} 

// Copy constructor
Index::Index(const Index& so) : sizes(so.sizes) {

	coord = new int[get_ndim()] ;
	for (int i=0 ; i<get_ndim() ; i++)
	   coord[i] = so.coord[i] ;
} 

// Constructor fro a tensor :
Index::Index (const Tensor& t) : sizes(t.valence) {
	for (int i=0 ; i<get_ndim() ; i++)
		sizes.set(i) = t.get_ndim() ;
	coord = new int[get_ndim()] ;
	for (int i=0 ; i<get_ndim() ; i++)
	    coord[i] = 0 ;
}

// Destructor
Index::~Index() {
	if(coord) delete [] coord ;
}

#ifdef ENABLE_MOVE_SEMANTIC
    Index::Index(Index && so) : sizes{std::move(so.sizes)}, coord{nullptr}
    {
        std::swap(coord,so.coord);
    }
    Index& Index::operator=(Index&& so)
    {
        sizes = std::move(so.sizes);
        std::swap(coord,so.coord);
        return *this;
    }
#endif

// Read/write
int& Index::set(int i) {
	assert(i>=0) ;
	assert(i<get_ndim()) ;
	return coord[i] ;
}

// Read only
int Index::operator() (int i) const {
	assert(i>=0) ;
	assert(i<get_ndim()) ;
	return coord[i] ;
}

// Assignement
void Index::operator= (const Index& so) {
	assert (sizes == so.sizes) ;
	for (int i=0 ; i<get_ndim() ; i++)
	    coord[i] = so.coord[i] ;
}

void Index::set_start() {
	for (int i=0 ; i<get_ndim() ; i++)
	     coord[i] = 0 ;
}

// Increment
bool Index::inc (int increm, int var) {
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

bool Index::operator== (const Index& xx) const {
	bool res = (get_ndim()==xx.get_ndim()) ? true : false ;
	if (res) 
		for (int i=0 ; i<get_ndim() ; i++)
			if (xx.coord[i] != coord[i])
				res = false ;
	return res ;
}

ostream& operator<< (ostream& o, const Index& so) {

	o << "(" ;
	for (int i=0 ; i<so.get_ndim()-1 ; i++)
	    o << so.coord[i] << ", " ;
	o << so.coord[so.get_ndim()-1] << ") in an array of " ;
	o << so.sizes << " points." ;
        return o ;
}
}
