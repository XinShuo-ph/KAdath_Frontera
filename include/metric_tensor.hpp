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

#ifndef __METRIC_TENSOR_HPP_
#define __METRIC_TENSOR_HPP_

#include "tensor.hpp"
#include "gsl/gsl_permutation.h"
#include <vector>
#include <cmath>

namespace Kadath {
double sign(gsl_permutation* p);
vector<int> ind_com(int i, int j, int k, int l);


/**
 * Particular type of \c Tensor, dedicated to the desription of metrics. It deals with symmetric, second order tensors only.
 * The two indices must be of the same type but can be either COV or CON.
 * \ingroup fields
 */
class Metric_tensor : public Tensor {

	public:
		/**
		* Constructor 
		* @param sp : the \c Space.
		* @param type_descr : the type of the indices (COV or CON)
		* @param ba : the tensorial basis.
		*/			
		Metric_tensor (const Space& sp, int type_descr, const Base_tensor&) ;
		/**
		* Constructor by copy.
		* If copie is false, the properties of the tensor are copied but not its values.
		*/
		Metric_tensor (const Metric_tensor&, bool copie = true) ; 
		Metric_tensor (const Space& sp, FILE*) ; ///< Copy from file.

		virtual ~Metric_tensor() ;

		// Assignement
        Metric_tensor & operator= (const Metric_tensor&) ; ///< Assignment to another \c Metric_tensor
		Metric_tensor & operator= (const Tensor& a) override ;
		Metric_tensor & operator= (double xx) override ;
 
		/**
		* Computes the inverse of the current objetc.
		* @return the inverse.
		*/
		Metric_tensor inverse();
	
		/**
		* Specialized base setting (for AADS)
		*/
		void std_base2(); 

		/**
		* @return the type of description (CON or COV)
		*/
		int get_type() const {return type_indice(0) ;} ;
} ;
}

#endif
