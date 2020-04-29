//
// Created by sauliac on 14/04/2020.
//
#include "config.h"

#include "memory.hpp"


namespace Kadath {
    Memory_mapper::mem_map_t Memory_mapper::memory_map{};
    Memory_mapper::ptr_list_t Memory_mapper::ptr_list{};
#ifdef HAVE_BOOST
    boost::pool<> Memory_mapper::size_4_memory_pool{4};
#endif
}

