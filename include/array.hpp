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


#include "headcpp.hpp"
#include "utilities.hpp"
#include "dim_array.hpp"
#include "index.hpp"
#include "array_iterator.hpp"

namespace Kadath {

/*!
 * Enum values to select how to write or read arrays from files.
 */
enum Array_ordering : bool {
    /*!
     * Value for first index major order (equivalent to column-major order for matrices).
     */
    first_index,
    /*!
     * Value for last index major order (raw-major order for matrices)
     */
    last_index
};

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
template <typename T> class Array : public Memory_mapped {
        static_assert(std::is_arithmetic<T>::value,"Array<T> implementation won't be efficient if T is not a "
                                                    " either an integral or floating point built-in type.");
	public:
	    Dim_array dimensions ;///< Dimensions of the \c Array.
	    int nbr ; ///< Total number of elements.
	    Memory_mapped_array<T> data ; ///< Elements of the \c Array.
	    
	public:
	//! Sylvain's stuff
        using reference = T& ;
        //! Sylvain's stuff
        using const_reference = T const &;
        //! Sylvain's stuff        
        using pointer = T*;
        //! Sylvain's stuff        
        using const_pointer = T const *;
	    /**
	    * Constructor from a \c Dim_array.
	    * The elements are not initialized.
	    */
	    explicit Array (const Dim_array& res) : dimensions{res}, nbr{1}, data{nullptr} {
            for (int i=0 ; i<res.get_ndim() ; i++) nbr *= res(i) ;
            data.resize(nbr) ;
        }
	     /**
	    * Constructor for a 1d-array.
	    * @param i [input] : size of the only dimension.
	    * The elements are not initialized.
	    */
	    explicit Array (int i)  : dimensions{1}, nbr{i}, data{nbr} {
             dimensions.set(0) = i ;
         }
   	    /**
	    * Constructor for a 2d-array.
	    * @param i [input] : size of the first dimension.
	    * @param j [input] : size of the second dimension.
	    * The elements are not initialized.
	    */
        explicit Array (int i, int j)  : dimensions{2}, nbr{i*j}, data{nbr} {
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
        explicit Array (int i, int j, int k)  : dimensions{3},nbr{i*j*k},data{nbr} {
            dimensions.set(0) = i ;
            dimensions.set(1) = j ;
            dimensions.set(2) = k ;
        }
	    /**
	     * Constructor from a file.
	     * The file should have been generated by the save function.
	     * @param  fd the file. 
	     * @param order specifies the type of index ordering used to write the array about to be read.
	     */
	    explicit Array (FILE* fd,Array_ordering order = last_index) ;

	public:
        //! Swaps contents between the two arrays (beware when using it with arrays of allocated pointers).
        void swap(Array<T> & so) {dimensions.swap(so.dimensions); std::swap(nbr,so.nbr); data.swap(so.data);}

	    void delete_data() {data.clear();} ///< Logical destructor (kills the data)

	    /**
	     * Resize the array by reallocating its ressource. All values are invalidated.
	     * @param new_dim new \c Dim_array object enumerating dimensions.
	     */
	    void resize(Dim_array const & new_dim) {this->operator=(std::move(Array<T>{new_dim}));}
	    /**
	     * Resize overload for the one-dimension array case.
	     * @param new_size new size of the array.
	     */
	    void resize(int new_size) {this->operator=(std::move(Array<T>{new_size}));}

	/**
	* Save in a file.
	* The file can then been used by the constructor from a file.
	* @param  fd the file. 
	* @param order allows to save the array in column major order (faster and need less memory, use the
	* \c first_index value - not compatible with files made using previous versions of Kadath) or row
	* major order (use the default \c last_index value). The most expensive case is used as default
	* for the sake of backward compatibility.
	*/
	    void save (FILE* fd, Array_ordering order = last_index) const ;
	    /**
	    * Assigns the same value to all the elements.
	    * @param xx [input] : value to be assigned.
	    */
	    Array & operator= (T xx) {for (int i=0 ; i<nbr ; i++) data[i] = xx ; return *this;}
	    /**
	    * Read/write of an element.
	    * @param pos [input] : position of the element.
	    */
        reference set(const Index& pos) {
	        assert (pos.sizes == dimensions) ;
            int index = pos(dimensions.get_ndim()-1) ;
            for (int i=dimensions.get_ndim()-2 ; i>=0 ; i--) {
                index *= dimensions(i) ;
		        assert ((pos(i) >=0) && (pos(i)<dimensions(i))) ;
                index += pos(i) ;
            }
            return data[index] ;
        }
//	    reference set(const Index& pos) {
////	assert (pos.sizes == dimensions) ;
//            int index = pos(0) ;
//            for (int i=1 ; i<dimensions.get_ndim() ; i++) {
//                index *= dimensions(i) ;
////		assert ((pos(i) >=0) && (pos(i)<dimensions(i))) ;
//                index += pos(i) ;
//            }
//            return data[index] ;
//        }
	   /**
	    * Read/write of an element.
	    * @param pos [input] : position of the element.
	    */
        reference set(const Array_iterator &pos) {return data[pos.position];}
	    /**
	    * Read/write of an element for a 1d-array.
	    * @param i [input] : position of the element.
	    */
	    reference set(int i) {
            assert (dimensions.get_ndim() == 1) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            return data[i] ;
        }
	    /**
	    * Read/write of an element for a 2d-array.
	    * @param i [input] : first index of the element.
	    * @param j [input] : second index of the element.
	    */
	    reference set(int i, int j)  {
            assert (dimensions.get_ndim() == 2) ;
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
            assert (dimensions.get_ndim() == 3) ;
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
            int index = pos(dimensions.get_ndim()-1) ;
            for (int i=dimensions.get_ndim()-2 ; i>=0 ; i--) {
                index *= dimensions(i) ;
		        assert ((pos(i) >=0) && (pos(i)<dimensions(i))) ;
                index += pos(i) ;
            }
            return data[index] ;
        }
//	    T operator() (const Index& pos) const {
////	assert (pos.sizes == dimensions) ;
//            int index = pos(0) ;
//            for (int i=1 ; i<dimensions.get_ndim() ; i++) {
//                index *= dimensions(i) ;
////		assert ((pos(i) >=0) && (pos(i)<dimensions(i))) ;
//                index += pos(i) ;
//            }
//            return data[index] ;
//	    }
 	    /**
	    * Read only of an element.
	    * @param pos [input] : position of the element.
	    */
	    T operator()(Array_iterator const &pos) const {return data[pos.position];}
	    /**
	    * Read only of an element for a 1d-array.
	    * @param i [input] : position of the element.
	    */
	    T operator() (int i) const {
            assert (dimensions.get_ndim() ==1) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            return data[i] ;
        }
	     /**
	    * Read only of an element for a 2d-array.
	    * @param i [input] : first index of the element.
	    * @param j [input] : second index of the element.
	    */
	    T operator() (int i,int j) const {
            assert (dimensions.get_ndim() ==2) ;
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
            assert (dimensions.get_ndim() ==3) ;
            assert ((i >=0) && (i<dimensions(0))) ;
            assert ((j >=0) && (j<dimensions(1))) ;
            assert ((k >=0) && (k<dimensions(2))) ;
            return data[i+dimensions(0)*(j+k*dimensions(1))] ;
//            return data[i*dimensions(1)*dimensions(2)+j*dimensions(2)+k] ;
        }
	    /**
	    * Direct accessor to the data, read only version
	    */
         CXX_17_ATTRIBUTES(nodiscard) const_pointer get_data() const {return data.get_data() ;} ;
	
	   /**
	   * Direct accessor to the data, read/write version
	   */
       CXX_17_ATTRIBUTES(nodiscard) pointer set_data() {return data.set_data() ;} ;

	    /**
	    * Returns the number of dimensions.
	    */
        CXX_17_ATTRIBUTES(nodiscard) int get_ndim() const {return dimensions.get_ndim() ;} ;
	    /**
	    * Returns the total number of elements.
	    */
        CXX_17_ATTRIBUTES(nodiscard) int get_nbr() const {return nbr ;} ;
	/**
	* Returns the size of a given dimension.
	*/
        CXX_17_ATTRIBUTES(nodiscard) int get_size(int i) const {return dimensions(i) ;} ;
	    /**
            * Returns the \c Dim_array of the \c Array.
	    */
	    CXX_17_ATTRIBUTES(nodiscard) const Dim_array& get_dimensions() const {return dimensions ;} ;
	    
	    /**
	    * Checks if a 1D array is increasing.
	    */
	    bool is_increasing() const ;

 	    /// Operator +
	    Array<T> & operator+= (const Array<T>& so)
        {assert(nbr == so.nbr); for(int i{0};i<nbr;i++) data[i] += so.data[i]; return *this;}
        /// Operator -=
	    Array<T> & operator-= (const Array<T>& so)
        {assert(nbr == so.nbr); for(int i{0};i<nbr;i++) data[i] -= so.data[i]; return *this;}
        /// Operator *=
	    Array<T> & operator*= (const Array<T>& so)
        {assert(nbr == so.nbr); for(int i{0};i<nbr;i++) data[i] *= so.data[i]; return *this;}
        /// Operator /=
	    Array<T> & operator/= (const Array<T>& so)
        {assert(nbr == so.nbr); for(int i{0};i<nbr;i++) data[i] /= so.data[i]; return *this;}
        /// Operator +=
	    Array<T> & operator+= (const T xx) {for(int i{0};i<nbr;i++) data[i] += xx; return *this;}  
	/// Operator -=
	    Array<T> & operator-= (const T xx) {for(int i{0};i<nbr;i++) data[i] -= xx; return *this;} 
	/// Operator *= 
	    Array<T> & operator*= (const T xx) {for(int i{0};i<nbr;i++) data[i] *= xx; return *this;}
	/// Operator /=
	    Array<T> & operator/= (const T xx) {for(int i{0};i<nbr;i++) data[i] /= xx; return *this;}

	    /**
	     * Index formulae to change major dimension indexation from first dimension to the last.
	     * WARNING : this method is not optimized and should not be called in called in computational loops, but rather
	     * only to import data from files.
	     */
	    int to_last_dim_major_index(int i_first_dim_major) const;
	    /**
	     * Index formulae to change major dimension indexation from last dimension to first.
	     * WARNING : this method is not optimized and should not be called in called in computational loops, but rather
	     * only to import data from files.
	     */
	    int to_first_dim_major_index(int i_last_dim_major) const;
	    
	friend class Matrice ;
      
        /// Output operator
	friend  ostream& operator<< <> (ostream&, const Array<T>& ) ;
	/// sine operator
	friend  Array<T> sin<> (const Array<T>&) ;
	/// cosine operator
	friend  Array<T> cos<> (const Array<T>&) ;
	/// Hyperbolic sine operator
	friend  Array<T> sinh<> (const Array<T>&) ;
	/// Hyperbolic sine operator
	friend  Array<T> cosh<> (const Array<T>&) ;
	/// Unitary operator +
	friend  Array<T> operator+ <>(const Array<T>&) ;
	/// Unitary operator -
	friend  Array<T> operator- <> (const Array<T>&) ;
	/// Operator + (two arrays)
	friend  Array<T> operator+ <> (const Array<T>&, const Array<T>&) ;
	/// Operator + (one array  and one value)
	friend  Array<T> operator+ <>(const Array<T>&, T) ;
	/// Operator + (one value  and one array)
	friend  Array<T> operator+ <>(T, const Array<T>&) ;
	/// Operator - (two arrays)
	friend  Array<T> operator- <>(const Array<T>&, const Array<T>&) ;
	/// Operator - (one array  and one value)
	friend  Array<T> operator- <>(const Array<T>&, T) ;
	/// Operator - (one value  and one array)
	friend  Array<T> operator- <>(T, const Array<T>&) ;
	/// Operator * (two arrays)
	friend  Array<T> operator* <>(const Array<T>&, const Array<T>&) ;
	/// Operator * (one array  and one value)
	friend  Array<T> operator* <>(const Array<T>&, T) ;
	/// Operator * (one value  and one array)
	friend  Array<T> operator* <>(T, const Array<T>&) ;
	/// Operator / (two arrays)
	friend  Array<T> operator/ <>(const Array<T>&, const Array<T>&) ;
	/// Operator / (one array  and one value)
	friend  Array<T> operator/ <>(const Array<T>&, T) ;
	/// Operator / (one value  and one array)
	friend  Array<T> operator/ <>(T, const Array<T>&) ;
	/// Operator power (integer version)
	friend  Array<T> pow <>(const Array<T>&, int) ;
	/// Operator power (double version)
	friend  Array<T> pow <>(const Array<T>&, double) ;
	/// Operator square root
	friend  Array<T> sqrt <>(const Array<T>&) ;
	/// Operator exponential
	friend  Array<T> exp <>(const Array<T>&) ;
	/// Operator logarithm
	friend  Array<T> log <>(const Array<T>&) ;
	/// Operator hyperbolic argtan
	friend  Array<T> atanh <>(const Array<T>&) ;
	/// Operator absolute value
	friend  Array<T> fabs <>(const Array<T>&) ;
	/// Scalar product of two arrays
	friend T scal<>(const Array<T>&, const Array<T>&) ;	
	/// Maximal difference between two arrays
	friend T diffmax<>(const Array<T>&, const Array<T>&) ;
	/// Maximal value
	friend T max<> (const Array<T>&) ;
	/// Minimal value
	friend T min<> (const Array<T>&) ;
	/// Summation of all the elements
	friend T sum<> (const Array<T>&) ;
	/// Operator arctan
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

    template<typename T> int Array<T>::to_last_dim_major_index(int i_first_dim_major) const {
        int const dim {dimensions.get_ndim()};
        int k{0};
        Index ret{dimensions};
        int q=i_first_dim_major;
        while(q!=0 && k<dim) {
            auto q_r = div(q,dimensions(k));
            q = q_r.quot;
            ret.set(k) = q_r.rem;
            assert(q==0 || k<dim-1);
            k++;
        }
        k=0;
        int ip = ret(k);
        for(;k<dim-1;k++) ip = ret(k+1) + ip * dimensions(k+1);
        return ip;
    }

    template<typename T> int Array<T>::to_first_dim_major_index(int i_last_dim_major) const {
        int const dim {dimensions.get_ndim()};
        int k{dim-1};
        Index ret{dimensions};
        int q = i_last_dim_major;
        while(q != 0 && k>=0) {
            auto q_r = div(q,dimensions(k));
            q = q_r.quot;
            ret.set(k) = q_r.rem;
            assert(q==0 || k>0);
            k--;
        }
        k = dim-1;
        int ip = ret(k);
        for(;k>0;k--) ip = ret(k-1) + dimensions(k-1)*ip;
        return ip;
    }

    template <typename T> Array<T>::Array (FILE* fd,Array_ordering order) : dimensions{fd},nbr{0},data{nullptr} {
        fread_be(&nbr, sizeof(int), 1, fd) ;
        data.resize(nbr) ;
        if(order == first_index) fread_be(data.set_data(), sizeof(T), nbr, fd) ;
        else {
            Memory_mapped_array<T> lim_data{nbr};
            fread_be(lim_data.set_data(),sizeof(T),nbr,fd);
            for(int i=0;i<nbr;i++) data[i] = lim_data[to_last_dim_major_index(i)];
        }
    }


    template <typename T> void Array<T>::save (FILE* fd,Array_ordering order) const {
        dimensions.save(fd) ;
        fwrite_be(&nbr, sizeof(int), 1, fd) ;
        if(order == first_index) fwrite_be(data.get_data(), sizeof(T), nbr, fd) ;
        else {
            Memory_mapped_array<T> lim_data {nbr};
            for(int i=0;i<nbr;i++) lim_data[i] = data[to_first_dim_major_index(i)];
            fwrite_be(lim_data.get_data(),sizeof(T),nbr,fd);
        }
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
