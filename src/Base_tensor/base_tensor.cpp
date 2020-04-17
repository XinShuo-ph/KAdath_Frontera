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

#include "headcpp.hpp"
#include "base_tensor.hpp"

namespace Kadath {
// Constructors

ostream& operator<< (ostream& flux , const Base_tensor& so) {
  flux << "Tensorial basis" << endl ;
  for (int d=0 ; d<so.space.get_nbr_domains() ; d++) {
    flux << "Domain " << d << " : " ;
    switch (so.get_basis(d)) {
      case -1 :
	flux << "Undefined" << endl ;
	break ;
      case CARTESIAN_BASIS :
	  flux << "Cartesian" << endl ;
	  break ;
      case SPHERICAL_BASIS :
	  flux << "Spherical" << endl ;
	  break ; 
	case MTZ_BASIS :
	  flux << "MTZ" << endl ;
	  break ;
      default:
	  cerr << "Unknown tensorial basis" << endl ;
	  abort() ;
    }
  }
  return flux ;
}

bool operator== (const Base_tensor& b1, const Base_tensor& b2) {
  assert(&b1.space==&b2.space) ;
  int res = true ;
  for (int i=0 ; i<b1.space.get_nbr_domains() ; i++)
    if (b1.get_basis(i) != b2.get_basis(i))
      res = false ;
  return res ;
}

}

