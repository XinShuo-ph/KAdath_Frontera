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

#include <cassert>
#include "headcpp.hpp"
#include "utilities.hpp"
#include "dim_array.hpp"
#include "index.hpp"
#include "array_iterator.hpp"

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
	    explicit Array (const Dim_array& res) : dimensions{res}, nbr{1}, data{nullptr} {
            for (int i=0 ; i<res.ndim ; i++) nbr *= res(i) ;
            data = new T[nbr] ;
        }
	     /**
	    * Constructor for a 1d-array.
	    * @param i [input] : size of the only dimension.
	    * The elements are not initialized.
	    */
	    explicit Array (int i)  : dimensions{1}, nbr{i}, data{new T[nbr]} {
             dimensions.set(0) = i ;
         }
   	    /**
	    * Constructor for a 2d-array.
	    * @param i [input] : size of the first dimension.
	    * @param j [input] : size of the second dimension.
	    * The elements are not initialized.
	    */
        explicit Array (int i, int j)  : dimensions{2}, nbr{i*j}, data{new T[nbr]} {
            dimensions.set(0) = i ;
            dimensions.set(1) = j ;
        }

	    /**
	    * Constructor for a 3d-array.
	    * @param i [input] : size of the first dimension.
	    * @param j [input] : size of the second dimension.
	    * @param k [input] : size of the third dimension.
	    * The elements are not initialized.
	    */
        explicit Array (int i, int j, int k)  : dimensions{3},nbr{i*j*k},data{new T[nbr]} {
            dimensions.set(0) = i ;
            dimensions.set(1) = j ;
            dimensions.set(2) = k ;
        }
	     /**
	    * Constructor from a file.
	    * The file should have been generated by the save function.
	    */
	    Array (FILE*) ;
	    //! Copy constructor.
	    Array (const Array<T>& so): dimensions{so.dimensions}, nbr{so.nbr}, data{new T[nbr]}
        {
            for (int i=0 ; i<nbr ; i++)
                data[i] = so.data[i] ;
        }
#ifdef ENABLE_MOVE_SEMANTIC
        //! Move constructor.
        Array (Array<T> && so) : dimensions{std::move(so.dimensions)}, nbr{so.nbr},
                              data{so.data} {so.data = nullptr;}
#endif
        //! Destructor.
	    ~Array() {delete_data();}

	public:
	  
	    void delete_data() {if(data) {delete [] data; data = nullptr;}} ///< Logical destructor (kills the data)
	  
		/**
	    * Save in a file.
	    * The file can then been used by the constructor from a file.
	    */
	    void save (FILE*) const ; 
	    //! Assignment operator.
	    void operator= (const Array<T>& so){
            for (int i=0 ; i<nbr ; i++)
            data[i] = so.data[i] ;
        }
#ifdef ENABLE_MOVE_SEMANTIC
        //! Move assignment operator.
        Array & operator= (Array<T> && so)
        {
            dimensions = std::move(so.dimensions);
            nbr = so.nbr;
            std::swap(data,so.data);
            return *this;
        }
#endif
	    /**
	    * Assigns the same value to all the elements.
	    * @param xx [input] : value to be assigned.
	    */
	    void operator= (T xx) {for (int i=0 ; i<nbr ; i++) data[i] = xx ;}
	    /**
	    * Read/write of an element.
	    * @param pos [input] : position of the element.
	    */

        reference set(const Index& pos) {
	        assert (pos.sizes == dimensions) ;
            int index = pos(dimensions.ndim-1) ;
            for (int i=dimensions.ndim-2 ; i>=0 ; i--) {
                index *= dimensions(i) ;
		        assert ((pos(i) >=0) && (pos(i)<dimensions(i))) ;
                index += pos(i) ;
            }
            return data[index] ;
        }
//	    reference set(const Index& pos) {
////	assert (pos.sizes == dimensions) ;
//            int index = pos(0) ;
//            for (int i=1 ; i<dimensions.ndim ; i++) {
//                index *= dimensions(i) ;
////		assert ((pos(i) >=0) && (pos(i)<dimensions(i))) ;
//                index += pos(i) ;
//            }
//            return data[index] ;
//        }
        reference set(const Array_iterator &pos) {return data[pos.position];}
	    /**
	    * Read/write of an element for a 1d-array.
	    * @param i [input] : position of the element.
	    */
	    reference set(int i) {
            assert (dimensions.ndim == 1) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            return data[i] ;
        }
	    /**
	    * Read/write of an element for a 2d-array.
	    * @param i [input] : first index of the element.
	    * @param j [input] : second index of the element.
	    */
	    reference set(int i, int j)  {
            assert (dimensions.ndim == 2) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            assert ((j >=0) && (j<dimensions(1))) ;
            return data[i+j*dimensions(0)] ;
//            return data[i*dimensions(1)+j] ;
        }
	    /**
	    * Read/write of an element for a 3d-array.
	    * @param i [input] : first index of the element.
	    * @param j [input] : second index of the element. 
            * @param k [input] : third index of the element.
	    */
	    reference set(int i, int j, int k) {
            assert (dimensions.ndim == 3) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            assert ((j >=0) && (j<dimensions(1))) ;
            assert ((k >=0) && (k<dimensions(2))) ;
            return data[i+dimensions(0)*(j+k*dimensions(1))] ;
//            return data[(i*dimensions(1)+j)*dimensions(2)+k] ;
        }
	    /**
	    * Read only of an element.
	    * @param pos [input] : position of the element.
	    */
        T operator() (const Index& pos) const {
            assert (pos.sizes == dimensions) ;
            int index = pos(dimensions.ndim-1) ;
            for (int i=dimensions.ndim-2 ; i>=0 ; i--) {
                index *= dimensions(i) ;
		        assert ((pos(i) >=0) && (pos(i)<dimensions(i))) ;
                index += pos(i) ;
            }
            return data[index] ;
        }
//	    T operator() (const Index& pos) const {
////	assert (pos.sizes == dimensions) ;
//            int index = pos(0) ;
//            for (int i=1 ; i<dimensions.ndim ; i++) {
//                index *= dimensions(i) ;
////		assert ((pos(i) >=0) && (pos(i)<dimensions(i))) ;
//                index += pos(i) ;
//            }
//            return data[index] ;
//	    }
	    T operator()(Array_iterator const &pos) const {return data[pos.position];}
	    /**
	    * Read only of an element for a 1d-array.
	    * @param i [input] : position of the element.
	    */
	    T operator() (int i) const {
            assert (dimensions.ndim ==1) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            return data[i] ;
        }
	     /**
	    * Read only of an element for a 2d-array.
	    * @param i [input] : first index of the element.
	    * @param j [input] : second index of the element.
	    */
	    T operator() (int i,int j) const {
            assert (dimensions.ndim ==2) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            assert ((j >=0) && (j<dimensions(1))) ;
             return data[i+j*dimensions(0)] ;
//            return data[i*dimensions(1)+j] ;
        }
		/**
	    * Read only of an element for a 3d-array.
	    * @param i [input] : first index of the element.
	    * @param j [input] : second index of the element. 
            * @param k [input] : third index of the element.
	    */
	    T operator() (int i,int j, int k) const {
            assert (dimensions.ndim ==3) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            assert ((j >=0) && (j<dimensions(1))) ;
            assert ((k >=0) && (k<dimensions(2))) ;
            return data[i+dimensions(0)*(j+k*dimensions(1))] ;
//            return data[i*dimensions(1)*dimensions(2)+j*dimensions(2)+k] ;
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

	    Array<T> & operator+= (const Array<T>& so)///< Operator +
        {assert(nbr == so.nbr); for(int i{0};i<nbr;i++) data[i] += so.data[i]; return *this;}
	    Array<T> & operator-= (const Array<T>& so)///< Operator -=
        {assert(nbr == so.nbr); for(int i{0};i<nbr;i++) data[i] -= so.data[i]; return *this;}
	    Array<T> & operator*= (const Array<T>& so)///< Operator *=
        {assert(nbr == so.nbr); for(int i{0};i<nbr;i++) data[i] *= so.data[i]; return *this;}
	    Array<T> & operator/= (const Array<T>& so)///< Operator /=
        {assert(nbr == so.nbr); for(int i{0};i<nbr;i++) data[i] /= so.data[i]; return *this;}
	    Array<T> & operator+= (const T xx) {for(int i{0};i<nbr;i++) data[i] += xx; return *this;}  ///< Operator +=
	    Array<T> & operator-= (const T xx) {for(int i{0};i<nbr;i++) data[i] -= xx; return *this;}  ///< Operator -=
	    Array<T> & operator*= (const T xx) {for(int i{0};i<nbr;i++) data[i] *= xx; return *this;}  ///< Operator *=
	    Array<T> & operator/= (const T xx) {for(int i{0};i<nbr;i++) data[i] /= xx; return *this;}  ///< Operator /=
	    
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

template <typename T> bool Array<T>::is_increasing() const {
        assert(get_ndim()==1) ;
        for (int i = 1 ; i < nbr ; i++){
            if (data[i]-data[i-1] <= 0){
                return false ;
            }
        }
        return true ;
    }


    template <typename T> Array<T>::Array (FILE* fd) : dimensions{fd},nbr{},data{nullptr} {
        fread_be(&nbr, sizeof(int), 1, fd) ;
        data = new T[nbr] ;
        fread_be(data, sizeof(T), nbr, fd) ;
    }


    template <typename T> void Array<T>::save (FILE* fd) const {
        dimensions.save(fd) ;
        fwrite_be(&nbr, sizeof(int), 1, fd) ;
        fwrite_be(data, sizeof(T), nbr, fd) ;
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
