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

#include "headcpp.hpp"
#include "dim_array.hpp"
#include "memory.hpp"

namespace Kadath {

class Tensor ;

/**
* Class that gives the position inside a multi-dimensional \c Array.
* It can also be used to give the indexes of a given component of a \c Tensor.
*
* It simply consists of a list of integers.
* \ingroup util
**/


class Index : public MemoryMappable {
    protected:
	/**
	* Sizes of the associated \c Array.
	* When used with a \c Tensor, it is the dimension, for each tensorial index.
	*/
	Dim_array sizes ;
	int* coord = nullptr ; ///< Value of each index.


    public:
	/**Standard constructor.
	* All the positions are set to zero.
	* @param dim [input] Sizes in each dimensions.
	**/
  explicit Index (const Dim_array& dim) ;
	Index (const Index&) ; ///< Constructor by copy
	Index (const Tensor&) ; ///< Constructor for looping on components of a tensor
	~Index() ; ///<Destructor.

	/**
	* Read/write of the position in a given dimension.
	* @param i [input] dimension.
	*/
	int& set(int) ;
	/**
	* Read/write of the position in a given dimension.
	* @param i [input] dimension.
	*/
	int operator() (int i) const ;
	/**
	* Returns the number of dimensions.
	*/
	int get_ndim() const {return sizes.get_ndim() ;} ;
	/**
	* Returns all the dimensions
	*/
	Dim_array const& get_sizes() const {return sizes ;} ;

	void set_start() ; ///< Sets the position to zero in all dimensions
	/**
	* Increments the position of the \c Index.
	* If one reaches the last point of a dimension, then the next one is increased.
	* @param increm [input] value of the increment.
	* @param var [input] dimension to be incremented.
	* @return \c false if the result is outside the \c Array and \c true otherwise.
	*/
	bool inc (int increm, int var=0) ;

	/**
	* Increments the position of the \c Index by one.
	* If one reaches the last point of a dimension, then the next one is increased.
	* @return \c false if the result is outside the \c Array and \c true otherwise.
	*/
	bool inc () ;

	void operator= (const Index&) ; ///< Assignement to annother \c Index.

	bool operator== (const Index&) const ; ///< Comparison operator

	template <class> friend class Array ;
	friend ostream& operator<< (ostream&, const Index&) ;
} ;
}

#include "implementation/index.cpp"

#endif
