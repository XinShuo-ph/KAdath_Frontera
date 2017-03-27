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

#include "ope_eq.hpp"
namespace Kadath {
Ope_eq::Ope_eq(const System_of_eqs* zesys, int dd, int nn) : syst(zesys), dom(dd), n_ope(nn) {
	parts = new Ope_eq* [n_ope] ;
	for (int i=0 ; i<n_ope ; i++)
		parts[i] = 0x0 ;
}


Ope_eq::Ope_eq(const System_of_eqs* zesys, int dd) : syst(zesys), dom(dd) {
}

Ope_eq::Ope_eq (const Ope_eq&) {
	cerr << "Copy constructor not defined for Ope_eq" << endl ;
	abort() ;
}

Ope_eq::~Ope_eq() {
	for (int i=0 ; i<n_ope ; i++)
		if (parts[i]!=0x0)
			delete parts[i] ;
	delete [] parts ;
}


}
