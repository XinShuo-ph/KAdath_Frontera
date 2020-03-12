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

#include "base_spectral.hpp"
#include "array.hpp"
#include "array.hpp"
#include "array_math.hpp"

namespace Kadath {
void coef_1d (int, Array<double>&) ;

void Base_spectral::coef_dim (int dim, int nbr_coef, Array<double> *& inout) const {

	Dim_array res_out(inout->get_dimensions()) ;
	res_out.set(dim) = nbr_coef ;
	Array<double> res (res_out) ;

	int after = 1 ;
	for (int i=0 ; i<dim ; i++)
		after *= inout->get_size(i) ;

	int before = 1 ;
	for (int i=dim+1 ; i<ndim ; i++)
	    before *= inout->get_size(i) ;

	int nbr_conf = inout->get_size(dim) ;
	int nbr = (nbr_coef > nbr_conf) ? nbr_coef : nbr_conf ;

	Index index_base (bases_1d[dim]->get_dimensions()) ;

	Index demarre_conf (inout->get_dimensions()) ;
	Index demarre_coef (res_out) ;

	Index loop_before_conf (inout->get_dimensions()) ;
	Index loop_before_out (res_out) ;

	Index lit_in (inout->get_dimensions()) ;
	Index put_out (res_out) ;
	
	Array<double> tab_1d (nbr) ;

	// Loop on dimensions before 
	for (int i=0 ; i<before ; i++) {
	    
	     demarre_conf = loop_before_conf ;
	     demarre_coef = loop_before_out ;
	    // On get la base
	    
	    int base = (*bases_1d[dim])(index_base) ;
	    // Loop on dimensions after :
	    for (int j=0 ; j<after ; j++) {
	    
	    	lit_in = demarre_conf ;
	    	for (int k=0 ; k<nbr_conf ; k++) {
	        	tab_1d.set(k) = (*inout)(lit_in) ;
			lit_in.inc(after) ;
			}
			
		// Transformation
		coef_1d(base, tab_1d) ;
		
		// On range :
		put_out = demarre_coef ;
		for (int k=0 ; k<nbr_coef ; k++) {
		    res.set(put_out) = tab_1d(k) ;
		    put_out.inc(after) ;
		}
		demarre_conf.inc() ;
		demarre_coef.inc() ;
        	}
		index_base.inc() ;
		loop_before_conf.inc(1, dim+1) ;
		loop_before_out.inc(1, dim+1) ;
	}
	
	delete inout ;
	inout = new Array<double>(res) ;
}

Array<double> Base_spectral::coef (const Dim_array& in_coef, const Array<double>& coloc) const {
	// Recopie
	Array<double>* cf = new Array<double> (coloc) ;
	for (int d=ndim-1 ; d>=0 ; d--) 
	      coef_dim(d, in_coef(d), cf) ;
	Array<double> res (*cf) ;
	delete cf ;
	return res ;
}

}

