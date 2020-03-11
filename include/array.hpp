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

#ifndef __ARRAY_HPP_
#define __ARRAY_HPP_

#include "assert.h"
#include "headcpp.hpp"
#include "index.hpp"
#include "dim_array.hpp"
namespace Kadath {
template <typename T> class Array ; 
template <typename T> ostream& operator<< (ostream&, const Array<T>& ) ;
template <typename T> Array<T> sin(const Array<T>& ) ;
template <typename T> Array<T> cos(const Array<T>& ) ;
template <typename T> Array<T> sinh(const Array<T>& ) ;
template <typename T> Array<T> cosh(const Array<T>& ) ;
template <typename T>  Array<T> operator+ (const Array<T>&) ;
template <typename T>  Array<T> operator- (const Array<T>&) ;
template <typename T>  Array<T> operator+ (const Array<T>&, const Array<T>&) ;
template <typename T>  Array<T> operator+ (const Array<T>&, T) ;
template <typename T>  Array<T> operator+ (T, const Array<T>&) ;
template <typename T>  Array<T> operator- (const Array<T>&, const Array<T>&) ;
template <typename T>  Array<T> operator- (const Array<T>&, T) ;
template <typename T>  Array<T> operator- (T, const Array<T>&) ;
template <typename T>  Array<T> operator* (const Array<T>&, const Array<T>&) ;
template <typename T>  Array<T> operator* (const Array<T>&, T) ;
template <typename T>  Array<T> operator* (T, const Array<T>&) ;
template <typename T>  Array<T> operator/ (const Array<T>&, const Array<T>&) ;
template <typename T>  Array<T> operator/ (const Array<T>&, T) ;
template <typename T>  Array<T> operator/ (T, const Array<T>&) ;
template <typename T>  Array<T> pow (const Array<T>&, int) ;
template <typename T>  Array<T> pow (const Array<T>&, double) ;
template <typename T>  Array<T> sqrt (const Array<T>&) ;
template <typename T>  Array<T> exp (const Array<T>&) ;
template <typename T>  Array<T> log (const Array<T>&) ;
template <typename T>  Array<T> atanh (const Array<T>&) ;
template <typename T>  Array<T> fabs (const Array<T>&) ;
template <typename T>  T scal (const Array<T>&, const Array<T>&) ;
template <typename T> T diffmax (const Array<T>&, const Array<T>&) ;
template <typename T>  T max (const Array<T>&) ;
template <typename T>  T min (const Array<T>&) ;
template <typename T>  T sum (const Array<T>&) ;
template <typename T>  Array<T> atan (const Array<T>&) ;

/**
* Template class for arrays.
*
* It is designed mainly to handle multi-dimensional arrays of \c int and \c double.
* \ingroup util
*/
template <typename T> class Array {
        static_assert(std::is_arithmetic<T>::value,"Array<T> implementation won't be efficient if T is not a "
                                                    " either an integral or floating point built-in type.");
	public:
	    Dim_array dimensions ;///< Dimensions of the \c Array.
	    int nbr ; ///< Total number of elements.
	    T * data ; ///< Elements of the \c Array.
	    
	public:
        using reference = T&;
        using const_reference = T const &;
        using pointer = T*;
        using const_pointer = T const *;

	    /**
	    * Constructor from a \c Dim_array.
	    * The elements are not initialized.
	    */
	    explicit Array (const Dim_array&) ;
	     /**
	    * Constructor for a 1d-array.
	    * @param i [input] : size of the only dimension.
	    * The elements are not initialized.
	    */
	    explicit Array (int i) ;
   	    /**
	    * Constructor for a 2d-array.
	    * @param i [input] : size of the first dimension.
	    * @param j [input] : size of the second dimension.
	    * The elements are not initialized.
	    */
            explicit Array (int i, int j) ;

	/**
	    * Constructor for a 3d-array.
	    * @param i [input] : size of the first dimension.
	    * @param j [input] : size of the second dimension.
	    * @param k [input] : size of the third dimension.
	    * The elements are not initialized.
	    */
            explicit Array (int i, int j, int k) ;
	 /**
	    * Constructor from a file.
	    * The file should have been generated by the save function.
	    */
	    Array (FILE*) ;
	    Array (const Array<T>&) ; ///< Copy constructor.
#ifdef ENABLE_MOVE_SEMANTIC
        Array (Array<T> &&) ; ///< Move constructor.
#endif
	    ~Array() ;	///<  Destuctor.

	public:
	  
	    void delete_data() {if(data) {delete [] data; data = nullptr;}} ///< Logical destructor (kills the data)
	  
		/**
	    * Save in a file.
	    * The file can then been used by the constructor from a file.
	    */
	    void save (FILE*) const ; 
	    void operator= (const Array<T>&) ; ///< Assignement to another \c Array
#ifdef ENABLE_MOVE_SEMANTIC
        Array & operator= (Array<T> &&) ; ///< Move assignment (transfer operator).
#endif
	    /**
	    * Assigns the same value to all the elements.
	    * @param xx [input] : value to be assigned.
	    */
	    void operator= (T xx) ;
	    /**
	    * Read/write of an element.
	    * @param pos [input] : position of the element.
	    */
	    reference set(const Index& pos) ;
	    /**
	    * Read/write of an element for a 1d-array.
	    * @param i [input] : position of the element.
	    */
	    reference set(int i) {
            /*assert (dimensions.ndim == 1) ;
            assert ((i >=0) && (i<dimensions(0))) ;*/
            return data[i] ;
        }
	    /**
	    * Read/write of an element for a 2d-array.
	    * @param i [input] : first index of the element.
	    * @param j [input] : second index of the element.
	    */
	    reference set(int i, int j)  {
            /*assert (dimensions.ndim == 2) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            assert ((j >=0) && (j<dimensions(1))) ;*/
            return data[i*dimensions(1)+j] ;
        }
	    /**
	    * Read/write of an element for a 3d-array.
	    * @param i [input] : first index of the element.
	    * @param j [input] : second index of the element. 
            * @param k [input] : third index of the element.
	    */
	    reference set(int i, int j, int k) {
            /*assert (dimensions.ndim == 3) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            assert ((j >=0) && (j<dimensions(1))) ;
            assert ((k >=0) && (k<dimensions(2))) ;*/
            return data[i*dimensions(1)*dimensions(2)+j*dimensions(2)+k] ;
        }
	    /**
	    * Read only of an element.
	    * @param pos [input] : position of the element.
	    */
	    T operator() (const Index&) const ;
	    /**
	    * Read only of an element for a 1d-array.
	    * @param i [input] : position of the element.
	    */
	    T operator() (int i) const {
            /*assert (dimensions.ndim ==1) ;
            assert ((i >=0) && (i<dimensions(0))) ;*/
            return data[i] ;
        }
	     /**
	    * Read only of an element for a 2d-array.
	    * @param i [input] : first index of the element.
	    * @param j [input] : second index of the element.
	    */
	    T operator() (int i,int j) const {
            /*assert (dimensions.ndim ==2) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            assert ((j >=0) && (j<dimensions(1))) ;*/
            return data[i*dimensions(1)+j] ;
        }
		/**
	    * Read only of an element for a 3d-array.
	    * @param i [input] : first index of the element.
	    * @param j [input] : second index of the element. 
            * @param k [input] : third index of the element.
	    */
	    T operator() (int i,int j, int k) const {
            /*assert (dimensions.ndim ==3) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            assert ((j >=0) && (j<dimensions(1))) ;
            assert ((k >=0) && (k<dimensions(2))) ;*/
            return data[i*dimensions(1)*dimensions(2)+j*dimensions(2)+k] ;
        }
	     /**
	   * Direct accessor to the data, read only version
	   */
	    const_pointer get_data() const {return data ;} ;
	
	   /**
	   * Direct accessor to the data, read/write version
	   */
	    pointer set_data() {return data ;} ;

	    /**
	    * Returns the number of dimensions.
	    */
	    int get_ndim() const {return dimensions.get_ndim() ;} ;
	    /**
	    * Returns the total number of elements.
	    */
	    int get_nbr() const {return nbr ;} ;
	    /**
	    * Returns the size of a given dimension.
	    * @param i [input] : dimension.
	    */
	    int get_size(int i) const {return dimensions(i) ;} ;
	    /**
            * Returns the \c Dim_array of the \c Array.
	    */
	    const Dim_array& get_dimensions() const {return dimensions ;} ;
	    
	    /*
	    * Checks if a 1D array is increasing.
	    */
	    bool is_increasing() const ;

	    void operator+= (const Array<T>&) ; ///< Operator +=
	    void operator-= (const Array<T>&) ; ///< Operator -=
	    void operator*= (const Array<T>&) ; ///< Operator *=
	    void operator/= (const Array<T>&) ; ///< Operator /=
	    void operator+= (const T&) ; ///< Operator +=
	    void operator-= (const T&) ; ///< Operator -=
	    void operator*= (const T&) ; ///< Operator *=
	    void operator/= (const T&) ; ///< Operator /=
	    
	friend class Matrice ;

	friend  ostream& operator<< <> (ostream&, const Array<T>& ) ;
	friend  Array<T> sin<> (const Array<T>&) ;
	friend  Array<T> cos<> (const Array<T>&) ;
	friend  Array<T> sinh<> (const Array<T>&) ;
	friend  Array<T> cosh<> (const Array<T>&) ;
	friend  Array<T> operator+ <>(const Array<T>&) ;
	friend  Array<T> operator- <> (const Array<T>&) ;
	friend  Array<T> operator+ <> (const Array<T>&, const Array<T>&) ;
	friend  Array<T> operator+ <>(const Array<T>&, T) ;
	friend  Array<T> operator+ <>(T, const Array<T>&) ;
	friend  Array<T> operator- <>(const Array<T>&, const Array<T>&) ;
	friend  Array<T> operator- <>(const Array<T>&, T) ;
	friend  Array<T> operator- <>(T, const Array<T>&) ;
	friend  Array<T> operator* <>(const Array<T>&, const Array<T>&) ;
	friend  Array<T> operator* <>(const Array<T>&, T) ;
	friend  Array<T> operator* <>(T, const Array<T>&) ;
	friend  Array<T> operator/ <>(const Array<T>&, const Array<T>&) ;
	friend  Array<T> operator/ <>(const Array<T>&, T) ;
	friend  Array<T> operator/ <>(T, const Array<T>&) ;
	friend  Array<T> pow <>(const Array<T>&, int) ;
	friend  Array<T> pow <>(const Array<T>&, double) ;
	friend  Array<T> sqrt <>(const Array<T>&) ;
	friend  Array<T> exp <>(const Array<T>&) ;	
	friend  Array<T> log <>(const Array<T>&) ;
	friend  Array<T> atanh <>(const Array<T>&) ;
	friend  Array<T> cos <>(const Array<T>&) ;
	friend  Array<T> sin <>(const Array<T>&) ;
	friend  Array<T> fabs <>(const Array<T>&) ;
	friend T scal<>(const Array<T>&, const Array<T>&) ;	
	friend T diffmax<>(const Array<T>&, const Array<T>&) ;
	friend T max<> (const Array<T>&) ;
	friend T min<> (const Array<T>&) ;
	friend T sum<> (const Array<T>&) ;
	friend  Array<T> atan <>(const Array<T>&) ;	
} ;

}
#endif
