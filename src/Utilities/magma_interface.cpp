/*
    Copyright 2020 sauliac

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

#include "magma_interface.hpp"

namespace Kadath{

#ifdef ENABLE_GPU_USE

    Magma_array::Magma_array(const Array<double> &source) : Base(source.get_nbr()),dim{static_cast<magma_int_t>(source.get_nbr())}
    {
        for(std::size_t i{0};i<this->size();i++)
        {
            std::size_t const k{i};
            (*this)[k] = source.get_data()[k];
        }
    }

    Magma_array & Magma_array::operator=(const Array<double> &source)
    {
        for(std::size_t i{0};i<this->size();i++)
        {
            std::size_t const k{i};
            (*this)[k] = source.get_data()[k];
        }
        return *this;
    }

    Magma_array & Magma_matrix::solve(Kadath::Magma_array &second_member) const
    {
        pivot.reset(new std::vector<magma_int_t,Magma_allocator<magma_int_t>>(order));
        magma_int_t info;
        magma_dgesv( order, 1, this->data(), lda, pivot->data(), second_member.data(), second_member.get_dim(), &info );
        return second_member;
    }


#endif
}