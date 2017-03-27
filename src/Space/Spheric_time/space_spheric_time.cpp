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
#include "spheric_time.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
namespace Kadath {
Space_spheric_time::Space_spheric_time(int ttype, const Dim_array& res, const Array<double>& bounds, double ttmin, double ttmax, bool wc) : tmin(ttmin), tmax (ttmax), withcompact (wc)  {

    // Verif :
    assert (bounds.get_ndim()== 1) ;

    ndim = 3 ;
    
    nbr_domains = (withcompact) ? bounds.get_size(0)+1 : bounds.get_size(0) ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;
    // Nucleus
    domains[0] = new Domain_spheric_time_nucleus(0, ttype, tmin, tmax, bounds(0), res) ;
    if (!withcompact) {
      for (int i=1 ; i<nbr_domains ; i++)
	domains[i] = new Domain_spheric_time_shell(i, ttype, tmin, tmax, bounds(i-1), bounds(i), res) ;
    }
    else {
      for (int i=1 ; i<nbr_domains-1 ; i++)
	domains[i] = new Domain_spheric_time_shell(i, ttype, tmin, tmax, bounds(i-1), bounds(i), res) ;
      domains[nbr_domains-1] = new Domain_spheric_time_compact (nbr_domains-1, ttype, tmin, tmax, bounds(nbr_domains-2), res) ;
    }
}

Space_spheric_time::Space_spheric_time(FILE* fd) {
	fread_be (&nbr_domains, sizeof(int), 1, fd) ;
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;
	fread_be (&tmin, sizeof(double), 1, fd) ;	
	fread_be (&tmax, sizeof(double), 1, fd) ;
	int indcomp ;
	fread_be (&indcomp, sizeof(int), 1, fd) ;
	withcompact = (indcomp==1) ? true : false ;
	domains = new Domain* [nbr_domains] ;
	//nucleus :
	domains[0] = new Domain_spheric_time_nucleus(0, fd) ;
	
	if (!withcompact) {
	//Shells :
	for (int i=1 ; i<nbr_domains ; i++)
	  domains[i] = new Domain_spheric_time_shell(i, fd) ;
	}
	else {
	  //Shells :
	for (int i=1 ; i<nbr_domains-1 ; i++)
	  domains[i] = new Domain_spheric_time_shell(i, fd) ;
	 domains[nbr_domains-1] = new Domain_spheric_time_compact(nbr_domains-1, fd) ;
	}
}

Space_spheric_time::~Space_spheric_time() {
}

void Space_spheric_time::save (FILE* fd) const  {
	fwrite_be (&nbr_domains, sizeof(int), 1, fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;	
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	fwrite_be (&tmin, sizeof(double), 1, fd) ;	
	fwrite_be (&tmax, sizeof(double), 1, fd) ;
	int indcomp = (withcompact) ? 1 : 0 ;
	fwrite_be (&indcomp, sizeof(int), 1, fd) ;
	for (int i=0 ; i<nbr_domains ; i++)
		domains[i]->save(fd) ;
}}
