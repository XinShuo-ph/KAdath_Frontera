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

#include "math.h"
#include "array.hpp"
  
namespace Kadath {
template <typename T> void Array<T>::operator+= (const Array<T>& so) {
	*this = *this + so ;
}

template <typename T> void Array<T>::operator-= (const Array<T>& so)  {
	*this = *this - so ;
}

template <typename T> void Array<T>::operator*= (const Array<T>& so) {
	*this = *this * so ;
}

template <typename T> void Array<T>::operator/= (const Array<T>& so) {
	*this = *this / so ;
}

template <typename T> void Array<T>::operator+= (const T& xx) {
	*this = *this + xx ;
}

template <typename T> void Array<T>::operator-= (const T& xx) {
	*this = *this - xx ;
}

template <typename T> void Array<T>::operator*= (const T& xx) {
	*this = *this * xx ;
}

template <typename T> void Array<T>::operator/= (const T& xx) {
	*this = *this/xx ;
}
	    
template <typename T> Array<T> sin (const Array<T>& so) {

	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = sin(so.data[i]) ;

	return res ;
}

template <typename T> Array<T> cos (const Array<T>& so) {

	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = cos(so.data[i]) ;

	return res ;
}

template <typename T>  Array<T> operator+ (const Array<T>& so) {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = so.data[i] ;

	return res ;
}


template <typename T>  Array<T> operator- (const Array<T>& so) {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = -so.data[i] ;

	return res ;
}


template <typename T>  Array<T> operator+ (const Array<T>& a, const Array<T>& b) {
	assert (a.dimensions==b.dimensions) ; 
	Array<T> res(a.dimensions) ;
	for (int i=0 ; i<a.nbr ; i++)
	    res.data[i] = a.data[i] + b.data[i] ;

	return res ;
}
	

template <typename T>  Array<T> operator+ (const Array<T>& so , T xx) {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = so.data[i] + xx;

	return res ;
}

template <typename T>  Array<T> operator+ (T xx, const Array<T>& so)  {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = xx + so.data[i] ;

	return res ;
}



template <typename T>  Array<T> operator- (const Array<T>& a, const Array<T>& b) {
	assert (a.dimensions==b.dimensions) ; 
	Array<T> res(a.dimensions) ;
	for (int i=0 ; i<a.nbr ; i++)
	    res.data[i] = a.data[i] - b.data[i] ;

	return res ;
}


template <typename T>  Array<T> operator- (const Array<T>& so, T xx) {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = so.data[i] - xx;

	return res ;
}


template <typename T>  Array<T> operator- (T xx, const Array<T>& so) {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = xx - so.data[i] ;

	return res ;
}

template <typename T>  Array<T> operator* (const Array<T>& a, const Array<T>& b) {
	assert (a.dimensions==b.dimensions) ; 
	Array<T> res(a.dimensions) ;
	for (int i=0 ; i<a.nbr ; i++)
	    res.data[i] = a.data[i] * b.data[i] ;

	return res ;
}


template <typename T>  Array<T> operator* (const Array<T>& so, T xx) {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = so.data[i] * xx;

	return res ;
}

template <typename T>  Array<T> operator* (T xx, const Array<T>& so) {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = xx * so.data[i] ;

	return res ;
}

template <typename T>  Array<T> operator/ (const Array<T>& a, const Array<T>& b ){
	assert (a.dimensions==b.dimensions) ; 
	Array<T> res(a.dimensions) ;
	for (int i=0 ; i<a.nbr ; i++)
	    res.data[i] = a.data[i] / b.data[i] ;

	return res ;
}

template <typename T>  Array<T> operator/ (const Array<T>& so, T xx)  {	
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = so.data[i] / xx;

	return res ;
}

template <typename T>  Array<T> operator/ (T xx, const Array<T>& so)  {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = xx / so.data[i] ;

	return res;
}

template <typename T>  Array<T> pow (const Array<T>& so, int n) {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = pow(so.data[i],n) ;

	return res;
}

template <typename T>  Array<T> pow (const Array<T>& so, double nn)  {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = pow(so.data[i],nn) ;

	return res;
}

template <typename T>  Array<T> sqrt (const Array<T>& so)  {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = sqrt(so.data[i]) ;

	return res;
}

template <typename T>  Array<T> exp (const Array<T>& so)  {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = exp(so.data[i]) ;

	return res;
}

template <typename T>  Array<T> log (const Array<T>& so)  {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = log(so.data[i]) ;

	return res;
}

template <typename T>  Array<T> fabs (const Array<T>& so)  {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = fabs(so.data[i]) ;

	return res;
}

template <typename T> T scal (const Array<T>& a, const Array<T>& b) {
	T res = a.data[0] * b.data[0] ;
	for (int i=1 ; i<a.get_size(0) ; i++)
		res += a.data[i] * b.data[i] ;
	return res ;
}


template <typename T> T diffmax (const Array<T>& a, const Array<T>& b) {
	assert (a.dimensions==b.dimensions) ;
	T res = 0 ;
	T diff ;
	for (int i=0 ; i<a.nbr ; i++) {
		diff = fabs(a.data[i]-b.data[i]) ;
		if (diff> res) 
			res = diff ;
	}
	return res ;
}

template <typename T>  T max (const Array<T>& so)  {
	T res = 0 ;
	for (int i=0 ; i<so.nbr ; i++)
		if (so.data[i]>res)
			res = so.data[i] ;
	return res;
}

template <typename T>  Array<T> atan (const Array<T>& so)  {
	Array<T> res(so.dimensions) ;
	for (int i=0 ; i<so.nbr ; i++)
	    res.data[i] = atan(so.data[i]) ;

	return res;
}
}
