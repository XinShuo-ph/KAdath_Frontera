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

#include "utilities.hpp"
#include "dim_array.hpp"
namespace Kadath {

Dim_array::Dim_array (FILE* fd) {
	fread_be(&ndim, sizeof(int), 1, fd) ;
	nbr = new int[ndim] ;
	fread_be(nbr, sizeof(int), ndim, fd) ;
}


void Dim_array::save (FILE* fd) const  {
	fwrite_be(&ndim, sizeof(int), 1, fd) ;
	fwrite_be(nbr, sizeof(int), ndim, fd) ;
}

// Output
ostream& operator<< (ostream& o, const Dim_array& so) {
	o << "(" ;
	for (int i=0 ; i<so.ndim-1 ; i++)
	    o << so.nbr[i] << ", " ;
	o << so.nbr[so.ndim-1] << ")" ;
        return o ;
}

}
