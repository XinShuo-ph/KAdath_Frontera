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

#ifndef __ARRAY_IMPL_HPP_
#define __ARRAY_IMPL_HPP_

#include "assert.h"
#include "array.hpp"
#include "headcpp.hpp"
#include "utilities.hpp"
#include "dim_array.hpp"

namespace Kadath {
// Constructor from a Dim_array
template <typename T> Array<T>::Array (const Dim_array& res) : dimensions(res){

	nbr = 1 ;
	for (int i=0 ; i<res.ndim ; i++)
	    nbr *= res(i) ;

	data = MemoryMapper::get_memory<T>(nbr);
}

// Constructor for a 1d-array
template <typename T> Array<T>::Array (int i) : dimensions(1){
	dimensions.set(0) = i ;
	nbr = i ;
	data = MemoryMapper::get_memory<T>(nbr);
}

template <typename T> Array<T>::Array (int i, int j) : dimensions(2){
	dimensions.set(0) = i ;
	dimensions.set(1) = j ;
	nbr = i*j ;
	data = MemoryMapper::get_memory<T>(nbr);
}

template <typename T> Array<T>::Array (int i, int j, int k) : dimensions(3){
	dimensions.set(0) = i ;
	dimensions.set(1) = j ;
	dimensions.set(2) = k ;
	nbr = i*j*k ;
	data = MemoryMapper::get_memory<T>(nbr);
}

// Copy constructor
template <typename T> Array<T>::Array (const Array<T>& so) : dimensions(so.dimensions){
	nbr = so.nbr ;
	data = MemoryMapper::get_memory<T>(nbr);

	for (int i=0 ; i<nbr ; i++)
	    data[i] = so.data[i] ;
}

template <typename T> Array<T>::Array (FILE* fd) : dimensions(fd) {
	fread_be(&nbr, sizeof(int), 1, fd) ;
	data = MemoryMapper::get_memory<T>(nbr);

	fread_be(data, sizeof(T), nbr, fd) ;
}

// Destructor
template <typename T> Array<T>::~Array() {
	delete_data() ;
}

template <typename T> void Array<T>::save (FILE* fd) const {
	dimensions.save(fd) ;
	fwrite_be(&nbr, sizeof(int), 1, fd) ;
	fwrite_be(data, sizeof(T), nbr, fd) ;
}

// Assignement
template <typename T> void Array<T>::operator= (const Array<T>& so) {
	assert (dimensions == so.dimensions) ;
	for (int i=0 ; i<nbr ; i++)
	    data[i] = so.data[i] ;
}

// Assignement to a given value
template <typename T> void Array<T>::operator= (T xx) {
	for (int i=0 ; i<nbr ; i++)
	    data[i] = xx ;
}

// Read/write
template <typename T> T& Array<T>::set(const Index& point) {
	assert (point.sizes == dimensions) ;
	int index = point(0) ;
	for (int i=1 ; i<dimensions.ndim ; i++) {
	        index *= dimensions(i) ;
		assert ((point(i) >=0) && (point(i)<dimensions(i))) ;
		index += point(i) ;
		}
	return data[index] ;
}

// Read/write for a 1d-array
template <typename T> T& Array<T>::set(int i) {
	assert (dimensions.ndim == 1) ;
	assert ((i >=0) && (i<dimensions(0))) ;
	return data[i] ;
}

// Read/write for a 2d-array
template <typename T> T& Array<T>::set(int i, int j) {
	assert (dimensions.ndim == 2) ;
	assert ((i >=0) && (i<dimensions(0))) ;
	assert ((j >=0) && (j<dimensions(1))) ;
	return data[i*dimensions(1)+j] ;
}

// Read/write for a 3d-array
template <typename T> T& Array<T>::set(int i, int j, int k) {
	assert (dimensions.ndim == 3) ;
	assert ((i >=0) && (i<dimensions(0))) ;
	assert ((j >=0) && (j<dimensions(1))) ;
	assert ((k >=0) && (k<dimensions(2))) ;
	return data[i*dimensions(1)*dimensions(2)+j*dimensions(2)+k] ;
}

// Read only
template <typename T>  T Array<T>::operator() (const Index& point) const {
	assert (point.sizes == dimensions) ;
	int index = point(0) ;
	for (int i=1 ; i<dimensions.ndim ; i++) {
	        index *= dimensions(i) ;
		assert ((point(i) >=0) && (point(i)<dimensions(i))) ;
		index += point(i) ;
		}

	return data[index] ;
}

// Read only for a 1d-array.
template <typename T> T Array<T>::operator() (int i) const {
	assert (dimensions.ndim ==1) ;
	assert ((i >=0) && (i<dimensions(0))) ;

	return data[i] ;
}
// Read only for a 2d-array.
template <typename T> T Array<T>::operator() (int i, int j) const {
	assert (dimensions.ndim ==2) ;
	assert ((i >=0) && (i<dimensions(0))) ;
	assert ((j >=0) && (j<dimensions(1))) ;
	return data[i*dimensions(1)+j] ;
}

// Read only for a 3d-array.
template <typename T> T Array<T>::operator() (int i, int j, int k) const {
	assert (dimensions.ndim ==3) ;
	assert ((i >=0) && (i<dimensions(0))) ;
	assert ((j >=0) && (j<dimensions(1))) ;
	assert ((k >=0) && (k<dimensions(2))) ;
	return data[i*dimensions(1)*dimensions(2)+j*dimensions(2)+k] ;
}

template <typename T> void Array<T>::delete_data() {
  if (data!= 0x0)
	  MemoryMapper::release_memory<T>(data,nbr);
  data=0x0 ;
}

// Display
template <typename T> ostream& operator<< (ostream& o, const Array<T>& so) {
	int ndim = so.dimensions.get_ndim() ;
	o << "Array of " << ndim << " dimension(s)" << endl ;
	Index xx (so.get_dimensions()) ;
	switch (ndim) {
		case 1:
		    for (int i=0 ; i<so.dimensions(0) ; i++) {
		        xx.set(0) = i ;
		        o << so(xx) << " " ;
			}
	            break ;
		case 2:
		    for (int i=0 ; i<so.dimensions(1) ; i++) {
		    	o << i << " : " ;
			xx.set(1) = i ;
			for (int j=0 ; j<so.dimensions(0) ; j++) {
			   xx.set(0) = j ;
			   o << so(xx) << " " ;
			   }
			o << endl ;
			}
		    break ;
		case 3:
		   for (int i=0 ; i<so.dimensions(2) ; i++) {
		    	o << "      " << i << endl ;
			xx.set(2) = i ;
			for (int j=0 ; j<so.dimensions(1) ; j++) {
		    	     o << j << " : " ;
			     xx.set(1) = j ;
			     for (int k=0 ; k<so.dimensions(0) ; k++) {
			          xx.set(0) = k ;
			          o << so(xx) << " " ;
				  }
			      o << endl ;
			     }
			o << endl ;
			}
		    break ;
		default:
		    for (int i=0 ; i<so.nbr ; i++)
		        o << so.data[i] << " " ;
	            o << endl ;
		    break ;
		}
	return o ;
}

}
#endif
