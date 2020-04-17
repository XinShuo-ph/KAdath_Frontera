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

#include "tensor.hpp"
namespace Kadath {

int add_m_quant (const Param_tensor & aa, const Param_tensor & bb) {
    if(aa && bb) {
        int val_a = aa.m_quant_affected ? aa.m_quant : 0 ;
        int val_b = bb.m_quant_affected ? bb.m_quant : 0 ;
        return (val_a<val_b) ? val_a : val_b ;
    }
    else return 0 ;
}

int mult_m_quant (const Param_tensor & aa, const Param_tensor & bb) {
    if(!aa && !bb) return 0;
    else if (!aa)  return (bb.m_quant_affected) ? bb.m_quant : 0 ;
    else if (!bb)  return (aa.m_quant_affected) ? aa.m_quant : 0 ;
    else {
        int val_a = (aa.m_quant_affected) ? aa.m_quant : 0;
        int val_b = (bb.m_quant_affected) ? bb.m_quant : 0;
        return (val_a + val_b);
    }
}
}
