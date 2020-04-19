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

Dim_array::Dim_array (FILE* fd) : Data_type{} {
    int read_size{};
	fread_be(&read_size, sizeof(int), 1, fd) ;
	this->resize(read_size) ;
	fread_be(this->set_data(), sizeof(int), size, fd) ;
}


void Dim_array::save (FILE* fd) const  {
    int const size{static_cast<int>(this->get_size())};
	fwrite_be(&size, sizeof(int), 1, fd) ;
	fwrite_be(this->get_data(), sizeof(int), size, fd) ;
}

// Output
ostream& operator<< (ostream& o, const Dim_array& so) {
	o << "(" ;
	for (int i=0 ; i<so.get_ndim()-1 ; i++)
	    o << so(i) << ", " ;
	o << so.back() << ")" ;
        return o ;
}

}
