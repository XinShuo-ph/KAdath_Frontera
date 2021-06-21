//
// This file is part of Margherita, the light-weight EOS framework
//
//  Copyright (C) 2017, Elias Roland Most
//                      <emost@th.physik.uni-frankfurt.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
#ifndef _HH_MARGHERITA_EOS
#define _HH_MARGHERITA_EOS


#include "cold_pwpoly.hh"
#include "cold_pwpoly_implementation.hh"

#include "cold_table.hh"
#include "cold_table_implementation.hh"

namespace Kadath {
namespace Margherita {

template class Cold_Table_t<0>;
using Cold_Table = Cold_Table_t<0>;

}
}
#endif
