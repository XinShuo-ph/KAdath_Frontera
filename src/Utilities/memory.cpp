/*
    Copyright 2019 Ludwig Jens Papenfort

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

#include "memory.hpp"

namespace Kadath {

MemoryMapper::mem_map_t MemoryMapper::mem_map;
MemoryMapper::ptr_list_t MemoryMapper::ptr_list;

#ifdef KADATH_VECTORMAP
MemoryMapper::mem_sizes_t MemoryMapper::mem_sizes;
#endif

}