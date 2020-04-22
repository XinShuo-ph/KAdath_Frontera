/*
    Copyright 2018 Philippe Grandclement

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


#include "array.hpp"
#include "array.hpp"
#include "list_comp.hpp"

namespace Kadath {

// Standard onsructor
List_comp::List_comp(int nc, int val) : ncomp(nc), valence(val) {
	pcomp = new Array<int>* [ncomp] ;
	for (int i=0 ; i<ncomp ; i++)
		pcomp[i] = new Array<int>(valence) ;
} 

// Constructor by copy
List_comp::List_comp(const List_comp& so) :  ncomp(so.ncomp), valence(so.valence) {
	pcomp = new Array<int>* [ncomp] ;
	for (int i=0 ; i<ncomp ; i++)
		pcomp[i] = new Array<int>(*so.pcomp[i]) ;
}

// Destructor
List_comp::~List_comp() {
	for (int i=0 ; i<ncomp ; i++)
		delete pcomp[i] ;
	delete [] pcomp ;
}

// Read/write
Array<int>* List_comp::set(int i) {
	assert(i>=0) ;
	assert(i<ncomp) ;
	return pcomp[i] ;
}


// Read only
Array<int>* List_comp::operator()(int i) const {
	assert(i>=0) ;
	assert(i<ncomp) ;
	return pcomp[i] ;
}
}

