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

#ifndef __DIM_ARRAY_HPP_
#define __DIM_ARRAY_HPP_

#include "headcpp.hpp"

namespace Kadath {
/**
* Class for storing the dimensions of an array
*
* It simply consists of a list of integers, being the get_size of a given \c Array, in each
* dimension.
* \ingroup util
**/

    class Dim_array {
    protected:
        int ndim ; ///< Number of dimensions.
        int* nbr ; ///< Size of the \c Array in each dimension.

    public:
        /** Standard constructor
        * @param nd [input] number of dimensions. The sizes are not initialized.
        **/
        explicit inline Dim_array (int dim) : ndim{dim}, nbr{new int[ndim]} {}
        inline Dim_array (const Dim_array &so): ndim{so.ndim}, nbr{new int[ndim]}
        {for (int i=0 ; i<ndim ; i++) nbr[i] = so.nbr[i] ;}
        ///< Constructor from a file (previously generated by the save member)
        inline ~Dim_array() {if(nbr) delete [] nbr;} ///< Destructor

        /**
        * Read/write of the get_size of a given dimension.
        * @param i [input] dimension.
        */
        inline int& set(int i) {assert(i>=0); assert(i<ndim); return nbr[i];}
        /**
        * Read only of the get_size of a given dimension.
        * @param i [input] dimension.
        */
        inline int operator() (int i) const {assert(i>=0); assert(i<ndim); return nbr[i];}
        /**
        * Returns the number of dimensions.
        */
        inline int get_ndim() const {return ndim ;} ;
        /**
         * Assignement to annother \c Dim_array.
         */
        inline void operator= (const Dim_array& so) {assert (ndim==so.ndim);for (int i=0 ;i<ndim;i++) nbr[i] = so.nbr[i];}

        inline void swap(Dim_array & so) noexcept {std::swap(ndim,so.ndim); std::swap(nbr,so.nbr);}
        ///< Save function

        template <class> friend class Array ;
        friend ostream& operator<< (ostream&, const Dim_array&) ;
        friend bool operator== (const Dim_array&, const Dim_array&) ;
        friend bool operator!= (const Dim_array&, const Dim_array&) ;

#ifdef ENABLE_MOVE_SEMANTIC
        inline Dim_array(Dim_array &&so) noexcept : ndim{so.ndim},nbr{so.nbr} {so.nbr = nullptr;} ///< Move constructor.
        inline Dim_array & operator=(Dim_array && so) noexcept {ndim = so.ndim; std::swap(nbr,so.nbr); return *this;}
#endif
    } ;

    ostream& operator<< (ostream&, const Dim_array&) ;
    inline bool operator== (const Dim_array& a, const Dim_array& b) {
        bool res = (a.ndim==b.ndim) ? true : false ;
        if (res)
            for (int i=0 ; i<a.ndim && res ; i++)
                res = (a.nbr[i] == b.nbr[i]);
        return res ;
    }
// Anti-comparison operator
    inline  bool operator!= (const Dim_array& a, const Dim_array& b) {return !(a==b) ;}
}
#endif
