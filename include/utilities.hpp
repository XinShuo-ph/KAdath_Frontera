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

#ifndef	__UTILITIES_HPP_
#define	__UTILITIES_HPP_



 
#include "stdio.h"
#include <string.h>
namespace Kadath {
class Param ;


/** Finding the zero a function.
 * 
 *  This routine locates a zero by means of the secant method. 
 * 
 *  @param (*f)(double, const Param\&) [input] Function the zero of which is 
 *		    to be searched: the routine computes x0 in a given
 *		    interval [a, b] such that 
 *			f(x0, par) = 0 , where par are the parameters of the
 *		    function f, stored in an object of the Lorene class 
 *		    \c Param . 
 *  @param par [input] Parameters of the function f.
 *  @param a [input] Lower bound of the search interval [a, b]
 *  @param b [input] Higher bound of the search interval [a, b]
 *  @param precis [input] Required precision in the determination of x0 : 
 *			the returned solution will be x0 +/- precis
 *  @param nitermax [input] Maximum number of iterations in the secant 
 *			    method to compute x0.
 *  @param niter [output] Number of iterations effectively used in computing x0				
 *  @return x0 (zero of function f)
 *
 */
double zerosec( double (*f)(double, const Param&), const Param& par, 
		double a, double b, double precis, int nitermax, 
		int& niter) ;


/** Writes integer(s) into a binary file according to the
 *  big endian convention.
 *
 *  This function has the same prototype and return value than
 *  the \c fwrite  function of the \c stdio  C library.
 *  The difference is that it ensures that the binary file is
 *  written in the big endian format, whatever the system is
 *  using little endian or big endian.
 *	@param aa [input] integer array to be written (in case of one
 *		element, address of this integer)
 *	@param size [input] number of bytes of one \c int  (must
 *		be 4)
 *	@param nb [input] number of elements in the array \c aa 
 *	@param fich [input] binary file (must have been
 *		open by \c fopen )
 *	@return number of integers effectively written in the file
 */		
int fwrite_be(const int* aa, int size, int nb, FILE* fich) ;

/** Writes double precision number(s) into a binary file according to the
 *  big endian convention.
 *
 *  This function has the same prototype and return value than
 *  the \c fwrite  function of the \c stdio  C library.
 *  The difference is that it ensures that the binary file is
 *  written in the big endian format, whatever the system is
 *  using little endian or big endian.
 *	@param aa [input] array of \c double  to be written (in case of one
 *		element, address of this \c double )
 *	@param size [input] number of bytes of one \c double  (must
 *		be 8)
 *	@param nb [input] number of elements in the array \c aa 
 *	@param fich [input] binary file (must have been
 *		open by \c fopen )
 *	@return number of \c double  effectively written in the file
 */		
int fwrite_be(const double* aa, int size, int nb, FILE* fich) ;

/** Reads integer(s) from a binary file according to the
 *  big endian convention.
 *
 *  This function has the same prototype and return value than
 *  the \c fread  function of the \c stdio  C library.
 *  The difference is that it assumes that the binary file is
 *  written in the big endian format, whatever the system is
 *  using little endian or big endian.
 *	@param aa [output] integer array to be read (in case of one
 *		element, address of this integer)
 *	@param size [input] number of bytes of one \c int  (must
 *		be 4)
 *	@param nb [input] number of elements in the array \c aa 
 *	@param fich [input] binary file (must have been
 *		open by \c fopen )
 *	@return number of integers effectively read in the file
 */		
int fread_be(int* aa, int size, int nb, FILE* fich) ;

/** Reads double precision number(s) from a binary file according to the
 *  big endian convention.
 *
 *  This function has the same prototype and return value than
 *  the \c fread  function of the \c stdio  C library.
 *  The difference is that it assumes that the binary file is
 *  written in the big endian format, whatever the system is
 *  using little endian or big endian.
 *	@param aa [output] array of \c double  to be read (in case of one
 *		element, address of this \c double )
 *	@param size [input] number of bytes of one \c double  (must
 *		be 8)
 *	@param nb [input] number of elements in the array \c aa 
 *	@param fich [input] binary file (must have been
 *		open by \c fopen )
 *	@return number of \c double  effectively read in the file
 */		
int fread_be(double* aa, int size, int nb, FILE* fich) ;
    
/** @} */


}
#endif
