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
#include "memory.hpp"

namespace Kadath {
/**
* Class for storing the dimensions of an array
*
* It simply consists of a list of integers, being the size of a given \c Array, in each
* dimension.
* \ingroup util
**/

class Dim_array  : public Memory_mapped_array<int> {
    public:
        using size_type = int;

        /** Standard constructor
        * @param nd [input] number of dimensions. The sizes are not initialized.
        **/
        explicit Dim_array (int dim) : Memory_mapped_array{dim} {}
        Dim_array (const Dim_array &so): Memory_mapped_array{so} {}
        Dim_array (FILE*) ; ///< Constructor from a file (previously generated by the save member)

        /**
        * Read/write of the size of a given dimension.
        * @param i [input] dimension.
        */
        int& set(int i) {assert(i>=0); assert(i<size); return data[i];}
        /**
        * Read only of the size of a given dimension.
        * @param i [input] dimension.
        */
        int operator() (int i) const {assert(i>=0); assert(i<size); return data[i];}
        /**
        * Returns the number of dimensions.
        */
        int get_ndim() const {return size ;} ;
        /**
         * Assignement to annother \c Dim_array.
         */
        void operator= (const Dim_array& so) {assert (size==so.size);for (int i=0 ;i<size;i++) data[i] = so.data[i];}

        void save (FILE*) const ; ///< Save function

#ifdef ENABLE_MOVE_SEMANTIC
        Dim_array(Dim_array &&so) : Memory_mapped_array<int>{std::forward<Dim_array&&>(so)} {}///< Move constructor.
        Dim_array & operator=(Dim_array && so) {Memory_mapped_array<int>::operator=(std::forward<Dim_array&&>(so)); return *this;}
#endif
} ;

ostream& operator<< (ostream&, const Dim_array&) ;
inline bool operator== (const Dim_array& a, const Dim_array& b) {
    bool res = (a.get_ndim()==b.get_ndim()) ? true : false ;
    if (res)
        for (int i=0 ; i<a.get_ndim() && res ; i++)
            res = (a(i) == b(i));
    return res ;
}
// Anti-comparison operator
inline  bool operator!= (const Dim_array& a, const Dim_array& b) {return !(a==b) ;}
}
#endif
