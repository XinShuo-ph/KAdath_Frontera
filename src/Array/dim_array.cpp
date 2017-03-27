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
#include "utilities.hpp"
#include "dim_array.hpp"
namespace Kadath {
// Standard onsructor
Dim_array::Dim_array(int dim) : ndim(dim) {
	nbr = new int[ndim] ;
} 

// Constructor by copy
Dim_array::Dim_array(const Dim_array& so) : ndim(so.ndim) {
	nbr = new int[ndim] ;
	for (int i=0 ; i<ndim ; i++)
	   nbr[i] = so.nbr[i] ;
}

Dim_array::Dim_array (FILE* fd) {
	fread_be(&ndim, sizeof(int), 1, fd) ;
	nbr = new int[ndim] ;
	fread_be(nbr, sizeof(int), ndim, fd) ;
}

// Destructor
Dim_array::~Dim_array() {
	delete [] nbr ;
}

// Read/write
int& Dim_array::set(int i) {
	assert(i>=0) ;
	assert(i<ndim) ;
	return nbr[i] ;
}

// Assignement
void Dim_array::operator= (const Dim_array& so) {
	assert (ndim==so.ndim) ;
	for (int i=0 ; i<ndim ; i++)
		nbr[i] = so.nbr[i] ;
}

// Read only
int Dim_array::operator() (int i) const {
	assert(i>=0) ;
	assert(i<ndim) ;
	return nbr[i] ;
}

void Dim_array::save (FILE* fd) const  {
	fwrite_be(&ndim, sizeof(int), 1, fd) ;
	fwrite_be(nbr, sizeof(int), ndim, fd) ;
}

// Output
ostream& operator<< (ostream& o, const Dim_array& so) {
	o << "(" ;
	for (int i=0 ; i<so.ndim-1 ; i++)
	    o << so.nbr[i] << ", " ;
	o << so.nbr[so.ndim-1] << ")" ;
        return o ;
}

// Comparison operator
bool operator== (const Dim_array& a, const Dim_array& b) {
	bool res = (a.ndim==b.ndim) ? true : false ;
	if (res)
	    for (int i=0 ; i<a.ndim ; i++)
	         if (a.nbr[i] != b.nbr[i])
		      res = false ;
	return res ;
}

// Anti-comparison operator
bool operator!= (const Dim_array& a, const Dim_array& b) {

	return !(a==b) ;
}
}
