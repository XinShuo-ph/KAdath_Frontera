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

#include "vector.hpp"
#include "tensor.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
namespace Kadath {
    Vector & Vector::operator=(const Vector& t) {
        assert (&espace==&t.espace) ;
        basis = t.basis ;
        assert(t.type_indice(0) == type_indice(0)) ;

        for (int i=0 ; i<3 ; i++) {
          *cmp[i] = *t.cmp[i] ;
        }
        return *this;
    }

    Vector & Vector::operator=(const Tensor& t) {

        assert (t.valence == 1) ;


        assert (&espace==&t.espace) ;
        basis = t.basis ;
        assert(t.type_indice(0) == type_indice(0)) ;

        for (int i=0 ; i<3 ; i++) {
          *cmp[i] = *t.cmp[i] ;
        }
        return *this;
    }

    Vector & Vector::operator=(double xx) {
        for (int i=0 ; i<3 ; i++) {
          *cmp[i] = xx ;
        }
        return *this;
    }

}