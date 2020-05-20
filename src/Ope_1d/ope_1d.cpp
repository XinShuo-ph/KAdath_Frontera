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
#include "val_domain.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"

namespace Kadath {
Array<double> Base_spectral::ope_1d (int (*func) (int, Array<double>&), 
			int var, const Array<double>& in, Base_spectral& base_out) const {
	
	Array<double> res (in.get_dimensions()) ;

	int after = 1 ;
	for (int i=0 ; i<var ; i++)
		after *= in.get_size(i) ;

	int before = 1 ;
	for (int i=var+1 ; i<ndim ; i++)
	    before *= in.get_size(i) ;

	int nbr = in.get_size(var) ;
	    
	Array_iterator index_base {*bases_1d[var]} ;

	Array_iterator demarre{in} ;
	Array_iterator loop_before{in} ;
	
	Array_iterator lit{in} ;
	Array_iterator put{in} ;
	
	Array<double> tab_1d (nbr) ;

	// Loop on dimensions before 
	for (int i=0 ; i<before ; i++) {
	    
	     demarre.set(loop_before) ;
	    // On get la base
	    
	    int base = (*bases_1d[var])(index_base) ;
	    // Loop on dimensions after :
	    for (int j=0 ; j<after ; j++) {
	    
	    	lit.set(demarre) ;
	    	for (int k=0 ; k<nbr ; k++) {
	        	tab_1d.set(k) = in(lit) ;
			    lit.inc(after) ;
            }
            // Transformation
            base_out.bases_1d[var]->set(index_base) = func(base, tab_1d) ;

            // On range :
            put.set(demarre) ;
            for (int k=0 ; k<nbr ; k++) {
                res.set(put) = tab_1d(k) ;
                put.inc(after) ;
            }

            demarre.inc() ;
        }
		index_base.inc() ;
		loop_before.inc1( var+1) ;
	}

	return res ;
}

}
