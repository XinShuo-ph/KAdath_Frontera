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
#include "index.hpp"
#include "tensor.hpp"

namespace Kadath {


// Constructor fro a tensor :
Index::Index (const Tensor& t) : sizes{t.valence}, coord{new int[get_ndim()]} {
	for (int i=0 ; i<get_ndim() ; i++) {
	    coord[i] = 0;
        sizes.set(i) = t.get_ndim();
    }
}

// Increment


ostream& operator<< (ostream& o, const Index& so) {

	o << "(" ;
	for (int i=0 ; i<so.get_ndim()-1 ; i++)
	    o << so.coord[i] << ", " ;
	o << so.coord[so.get_ndim()-1] << ") in an array of " ;
	o << so.sizes << " points." ;
        return o ;
}
}
