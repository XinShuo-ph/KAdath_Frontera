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

int Param_tensor::get_m_order() const {
  if (!m_order_affected) {
    cerr << "m_order not affected in Param_tensor" << endl ;
    abort() ;
  }
  else
    return m_order ;
}

int Param_tensor::get_m_quant() const {
  if (!m_quant_affected) {
    cerr << "m_quant not affected in Param_tensor" << endl ;
    abort() ;
  }
  else
    return m_quant ;
}

int& Param_tensor::set_m_quant() {
  m_quant_affected = true ;
  return m_quant ;
}

int& Param_tensor::set_m_order() {
  m_order_affected = true ;
  return m_order ;
}

int add_m_quant (const Param_tensor* aa, const Param_tensor* bb) {
  int res = 0 ;
  if ((aa==0x0) || (bb==0x0))
    res = 0 ;
  else {
    int val_a = (aa->m_quant_affected) ? aa->m_quant : 0 ;
    int val_b = (bb->m_quant_affected) ? bb->m_quant : 0 ;
    res = (val_a<val_b) ? val_a : val_b ;
  }
  return res ;
}

int mult_m_quant (const Param_tensor* aa, const Param_tensor* bb) {
  int res = 0 ;
  if ((aa==0x0) && (bb==0x0))
    return 0 ;
  if (aa==0x0) {
      res = (bb->m_quant_affected) ? bb->m_quant : 0 ;
      return res ;
  }
  if (bb==0x0)  {
      res = (aa->m_quant_affected) ? aa->m_quant : 0 ;
      return res ;
  }
  int val_a = (aa->m_quant_affected) ? aa->m_quant : 0 ;
  int val_b = (bb->m_quant_affected) ? bb->m_quant : 0 ;
  return (val_a+val_b) ;
}
}
