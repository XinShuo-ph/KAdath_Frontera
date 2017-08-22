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
#include "polar_periodic.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
namespace Kadath {
Space_polar_periodic::Space_polar_periodic(int ttype, double omega, const Dim_array& res, const Array<double>& bounds) {

    // Verif :
    assert (bounds.get_ndim()== 1) ;

    ndim = 3 ;
    
    nbr_domains = bounds.get_size(0) ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
    // Nucleus
    domains[0] = new Domain_polar_periodic_nucleus(0, ttype, bounds(0), omega, res) ;
    for (int i=1 ; i<nbr_domains ; i++)
       domains[i] = new Domain_polar_periodic_shell(i, ttype, bounds(i-1), bounds(i), omega, res) ;
  
}

Space_polar_periodic::Space_polar_periodic(FILE* fd) {
	fread_be (&nbr_domains, sizeof(int), 1, fd) ;
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;
	domains = new Domain* [nbr_domains] ;
	//nucleus :
	domains[0] = new Domain_polar_periodic_nucleus(0, fd) ;
	//Shells :
	for (int i=1 ; i<nbr_domains-1 ; i++)
		domains[i] = new Domain_polar_periodic_shell(i, fd) ;
}

Space_polar_periodic::~Space_polar_periodic() {
}

void Space_polar_periodic::save (FILE* fd) const  {
	fwrite_be (&nbr_domains, sizeof(int), 1, fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;	
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	for (int i=0 ; i<nbr_domains ; i++)
		domains[i]->save(fd) ;
}
}
