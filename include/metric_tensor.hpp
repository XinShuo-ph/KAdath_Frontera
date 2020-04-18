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
    int metric_tensor_position_array (const Array<int>& idx, int ndim);
    int metric_tensor_position_index (const Index& idx, int ndim);
    Array<int> metric_tensor_indices (int pos, int valence, int ndim);

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

            // Assignement
            Metric_tensor & operator= (const Metric_tensor&) ; ///< Assignment to another \c Metric_tensor
            Metric_tensor & operator= (const Tensor& a) override ;
            Metric_tensor & operator= (double xx) override ;

#ifdef TENSOR_MOVE_SEMANTIC
            Metric_tensor(Metric_tensor &&s) noexcept : Tensor{std::move(s)} {}
            Metric_tensor & operator=(Metric_tensor &&s) noexcept {this->Tensor::operator=(std::move(s)); return *this;}
#endif

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

    inline Metric_tensor::Metric_tensor (const Space& sp, int type_descr, const Base_tensor& bb) :
            Tensor (sp, 2, type_descr, 6, bb, 3) {

        give_place_array = metric_tensor_position_array ;
        give_place_index = metric_tensor_position_index ;
        give_indices = metric_tensor_indices ;
    }

    inline Metric_tensor::Metric_tensor (const Metric_tensor& so, bool copie) :
            Tensor (so, copie) {
    }

    inline Metric_tensor::Metric_tensor (const Space& sp, FILE* ff) :
            Tensor (sp, 3, ff) {

        assert (valence==2) ;
        assert (type_indice(0)==type_indice(1)) ;
        assert (n_comp==6) ;

        // Overwrite the storage functions :
        give_place_array = metric_tensor_position_array ;
        give_place_index = metric_tensor_position_index ;
        give_indices = metric_tensor_indices ;
    }

}

#endif
