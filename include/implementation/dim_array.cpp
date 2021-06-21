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

#ifndef __DIM_ARRAY_IMPL_CPP_
#define __DIM_ARRAY_IMPL_CPP_

#include "assert.h"
#include "utilities.hpp"

namespace Kadath {


// Standard onsructor
inline Dim_array::Dim_array(int dim) : ndim(dim) {
	nbr = MemoryMapper::get_memory<int>(ndim);
}

// Constructor by copy
inline Dim_array::Dim_array(const Dim_array& so) : ndim(so.ndim) {
	nbr = MemoryMapper::get_memory<int>(ndim);
	for (int i=0 ; i<ndim ; i++)
	   nbr[i] = so.nbr[i] ;
}

inline Dim_array::Dim_array (FILE* fd) {
	fread_be(&ndim, sizeof(int), 1, fd) ;
	nbr = MemoryMapper::get_memory<int>(ndim);
	fread_be(nbr, sizeof(int), ndim, fd) ;
}

// Destructor
inline Dim_array::~Dim_array() {
  MemoryMapper::release_memory<int>(nbr, ndim);
}

// Read/write
inline int& Dim_array::set(int i) {
	assert(i>=0) ;
	assert(i<ndim) ;
	return nbr[i] ;
}

// Assignement
inline void Dim_array::operator= (const Dim_array& so) {
	assert (ndim==so.ndim) ;
	for (int i=0 ; i<ndim ; i++)
		nbr[i] = so.nbr[i] ;
}

// Read only
inline int Dim_array::operator() (int i) const {
	assert(i>=0) ;
	assert(i<ndim) ;
	return nbr[i] ;
}

inline void Dim_array::save (FILE* fd) const  {
	fwrite_be(&ndim, sizeof(int), 1, fd) ;
	fwrite_be(nbr, sizeof(int), ndim, fd) ;
}
}

#endif