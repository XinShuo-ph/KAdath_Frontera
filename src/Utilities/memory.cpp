//
// Created by sauliac on 14/04/2020.
//
#include "config.h"

#include "memory.hpp"

namespace Kadath {
    namespace Memory {
        template<> Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::Data
                Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::data{};
    }
}

