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

#ifndef __POINT_HPP_
#define __POINT_HPP_

#include "headcpp.hpp"
namespace Kadath {
/**
 * The class \c Point is used to store the coordinates of a point.
 * \ingroup fields
 */

class Point {
    protected:
        int ndim ; ///< Number of dimensions.
	double* coord ; ///< Array on the coordinates (mainly designed for absolute Cartesian coordinates).

    public:
	/**
	* Standard constructor (the coordinates are not affected).
	* @param n [input] : number of dimensions.
	*/
        explicit Point (int n) ;
        Point (FILE*) ; ///< Constructor from a file
	Point (const Point&) ; ///< Constructor by copy.
	~Point() ; ///< Destuctor

#ifdef ENABLE_MOVE_SEMANTIC
    Point(Point &&);///< Move constructor.
    Point & operator=(Point &&);///<Move assignment.
#endif
	
	void save (FILE*) const ; ///< Saving function
	void operator= (const Point&) ; ///< Assignement to another \c Point
	double& set(int) ; ///< Read/write of a coordinate
	double operator() (int) const ; ///< Read only of a coordinate.
	/**
	* @returns the number of dimensions.
	*/
	const int& get_ndim() const {return ndim ;} ; 
	
	friend ostream& operator<< (ostream&, const Point&) ; ///< Display
} ; 
}
#endif
