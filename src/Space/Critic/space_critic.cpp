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
#include "critic.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
namespace Kadath {
Space_critic::Space_critic(int ttype, double xlim, const Dim_array& res_inner, const Dim_array& res_outer) {

    ndim = 2 ;
    nbr_domains = 2 ;
    type_base = ttype ;
    // Two domains
    domains = new Domain* [2] ;
    // Inner one
    domains[0] = new Domain_critic_inner(0, ttype, res_inner, xlim) ;
    // Outer one
    domains[1] = new Domain_critic_outer(1, ttype, res_outer, xlim) ;
}

Space_critic::Space_critic(FILE* fd) {
	nbr_domains = 2 ;
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;
	domains = new Domain* [2] ;
	// Inner
	domains[0] = new Domain_critic_inner(0, fd) ;
	// Outer
	domains[1] = new Domain_critic_outer(1, fd) ;
}

Space_critic::~Space_critic() {
}

void Space_critic::save (FILE* fd) const  {
	fwrite_be (&ndim, sizeof(int), 1, fd) ;	
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	domains[0]->save(fd) ;
	domains[1]->save(fd) ;
}

Array<int> Space_critic::get_indices_matching_non_std(int dom, int bound) const {

	switch (dom)  {
		case 0 : {
			// Inner
			Array<int> res (2,1) ;
			switch (bound) {
				case OUTER_BC :
					res.set(0,0) = 1 ; // Outer domain
					res.set(1,0) = INNER_BC ;
					break ;
				default :
					cerr << "Bad bound in Space_critic::get_indices_matching_non_std" << endl ;
					abort() ;
				}
			return res ;
			}
		case 1 : 
			{
			// Outer domain
			Array<int> res(2, 1) ;
			switch (bound) {
				case INNER_BC :
					res.set(0,0) = 0 ; // Inner domain
					res.set(1,0) = OUTER_BC ;
					break ;
				default :
					cerr << "Bad bound in Space_critic::get_indices_matching_non_std" << endl ;
					abort() ;
				}
			return res ;
			}
		default :
			cerr << "Bad domain in Space_critic::get_indices_matching_non_std" << endl ;
			abort() ;
		}
}
}
