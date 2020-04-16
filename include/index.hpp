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

#ifndef __INDEX_HPP_
#define __INDEX_HPP_

#include <cassert>
#include "headcpp.hpp"
#include "dim_array.hpp"

namespace Kadath {

class Tensor;

/**
* Class that gives the position inside a multi-dimensional \c Array.
* It can also be used to give the indexes of a given component of a \c Tensor.
*
* It simply consists of a list of integers.
* \ingroup util
**/
class Index MEMORY_ALLOCATION_POLICY {
protected:
	/**
	* Sizes of the associated \c Array.
	* When used with a \c Tensor, it is the dimension, for each tensorial index.
	*/
	Dim_array sizes ;
	Memory_mapped_array<int> coord ; ///< Value of each index.


public:
	/**Standard constructor.
	* All the positions are set to zero.
	* @param dim [input] Sizes in each dimensions.
	**/
    explicit Index (const Dim_array& dim) : sizes{dim}, coord{dim.get_ndim()}
    {for(unsigned i{0u};i<dim.get_ndim();i++) coord[i]=0;}
    /**
     * Copy constructor.
     * @param so source to be copied.
     */
	Index (const Tensor&) ; ///< Constructor for looping on components of a tensor

	void swap(Index & so) {sizes.swap(so.sizes); coord.swap(so.coord);}

	/**
	* Read/write of the position in a given dimension.
	* @param i [input] dimension.
	*/
	int& set(int i) {assert(i>=0); assert(i<get_ndim()); return coord[i] ;}
	/**
	* Read/write of the position in a given dimension.
	* @param i [input] dimension.
	*/
	int operator() (int i) const {assert(i>=0); assert(i<get_ndim()); return coord[i]; }
	/**
	* Returns the number of dimensions.
	*/
	int get_ndim() const {return sizes.get_ndim() ;} ;
	/**
	* Returns all the dimensions
	*/
	Dim_array const& get_sizes() const {return sizes ;} ;
	
	void set_start()  ///< Sets the position to zero in all dimensions
    {
        for (int i=0 ; i<get_ndim() ; i++) coord[i] = 0 ;
    }
	/**
	* Increments the position of the \c Index.
	* If one reaches the last point of a dimension, then the next one is increased.
	* @param increm [input] value of the increment.
	* @param var [input] dimension to be incremented.
	* @return \c false if the result is outside the \c Array and \c true otherwise.
	*/
	bool inc (int increm, int var=0) {
        int const ndimm1{get_ndim()-1};
        if((var >=0) && (var<=ndimm1)) {
            coord[var] += increm ;
            div_t division {div(coord[var],sizes(var))};
            for (int i=var ; i<ndimm1 && division.quot>0; i++) {
                division = div(coord[i], sizes(i)) ;
                coord[i] = division.rem ;
                coord[i+1] += division.quot ;
            }
            return (coord[ndimm1] < sizes(ndimm1));
        }
        else return false;
    }
    bool inc1(int var)  {
        int const ndimm1{get_ndim()-1};
        int i{var};
        while(i<= ndimm1 && coord[i]==(sizes(i)-1))
        {
            coord[i] = 0;
            i++;
        }
        if(i>=(ndimm1+1)) { return false;}
        else {
            coord[i]++;
            return true;
        }
	}
	bool inc() {return inc1(0);}


	// Increment
    bool inc_vanilla (int increm, int var) {
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

        /**
	 * Assignment operator.
	 * @param so source to copy from.
	 */
	void operator=(const Index& so) {assert(sizes==so.sizes);for(int i=0 ;i<get_ndim();i++)coord[i] = so.coord[i];}
	
	bool operator== (const Index& xx) const {
        bool res = (get_ndim()==xx.get_ndim()) ;
        if (res)
            for (int i=0 ; i<get_ndim() && res; i++)
                res = (xx.coord[i] == coord[i]);
        return res ;
    }

	template <class> friend class Array ;
	friend ostream& operator<< (ostream&, const Index&) ;
} ;
}
#endif
