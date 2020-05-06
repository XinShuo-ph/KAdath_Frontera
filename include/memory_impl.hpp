//
// Created by sauliac on 30/04/2020.
// This file is designed to be #included just after the declaration of the Memory_mapper template class in
// in memory.hpp. Iknow, this is ugly design but it keeps the header less verbose and the method inlinable.
// DO NOT ATTEMPT TO INCLUDE IT ANYWHERE ELSE.
//

template<Map_type_tag mtt> void * Memory_mapper<mtt,Memory_pools_use::disabled>::get_memory(std::size_t const sz) {
    if (sz == 0) return nullptr;
    else {
        void *raw_mem_ptr;
        //first find the entry of the size sz in the map (or create it if it doesn't exists)
#ifdef MAP_MEMORY_WITH_VECTOR
        auto pos = std::find_if(data.memory_map.begin(), data.memory_map.end(),
                                            [sz](std::pair<std::size_t, ptr_vec_t> const &x) { return x.first == sz; });
        if (pos == data.memory_map.end()) {
            data.memory_map.emplace_back(sz, ptr_vec_t{});
            pos = std::prev(data.memory_map.end());
        }
        auto &mem = pos->second;
//        if (pos == data.memory_map.end()) data.memory_map.emplace_back(sz,ptr_vec_t{});
//        std::size_t index = pos - mem_sizes.begin();
//        auto &mem = data.memory_map[index];
#else
        auto &mem = data.memory_map[sz];
#endif
        if (mem.empty()) {
            data.ptr_list.emplace_back(std::malloc(sz));
            raw_mem_ptr = data.ptr_list.back().get();
        } else {
            raw_mem_ptr = mem.back();
            mem.pop_back();
        }
        return raw_mem_ptr;
    }
}

template<Map_type_tag mtt> void Memory_mapper<mtt,Memory_pools_use::disabled>::release_memory(void *raw_mem_ptr, std::size_t const sz) {
    if (raw_mem_ptr != nullptr) {
#ifdef MAP_MEMORY_WITH_VECTOR
        auto pos = std::find_if(data.memory_map.begin(), data.memory_map.end(),
                                [sz](std::pair<std::size_t, ptr_vec_t> const &x) {
                                    return x.first == sz;
                                });
#ifdef ALL_CHECKS_ENABLED
        if (pos == data.memory_map.end()) {
            std::cerr << "Error : Memory_mapper::release_memory(p,s) with p = @" << raw_mem_ptr
                      << " and s=" << sz << ".\n Attempting to free unallocated memory chunk.\n";
            abort();
        }
#endif // ifdef ALL_CHECKS_ENABLED
        pos->second.push_back(raw_mem_ptr);
#else //ifdef MAP_MEMORY_WITH_VECTOR
#ifdef ALL_CHECKS_ENABLED
        auto pos = data.memory_map.find(sz);
        if (pos == data.memory_map.end()) {
            std::cerr << "Error : Memory_mapper::release_memory(p,s) with p = @" << raw_mem_ptr
                      << " and s=" << sz <<".\n Attempting to free unallocated memory chunk.\n";
            abort();
        }
        else {
            pos->second.push_back(raw_mem_ptr);
        }
#else  //ifdef ALL_CHECKS_ENABLED
        data.memory_map[sz].push_back(raw_mem_ptr);
#endif //ifdef ALL_CHECKS_ENABLED
#endif //ifdef MAP_MEMORY_WITH_VECTOR
    }
}

#ifdef HAVE_BOOST

template<Map_type_tag mtt> void * Memory_mapper<mtt,Memory_pools_use::enabled>::get_memory(std::size_t const sz) {
#ifdef MAP_MEMORY_WITH_VECTOR
    auto pos = std::find_if(data.memory_map.begin(), data.memory_map.end(),
                            [sz](std::pair<std::size_t, Memory_pool> const &x) { return x.first == sz; });
#else //ifdef MAP_MEMORY_WITH_VECTOR
    auto pos = data.memory_map.find(sz);
#endif //ifdef MAP_MEMORY_WITH_VECTOR
    if (pos == data.memory_map.end()) {
        data.memory_map.emplace_back(sz, Memory_pool{new boost::pool<>{sz}});
        pos = std::prev(data.memory_map.end());
    }
    return pos->second->malloc();
}

template<Map_type_tag mtt> void Memory_mapper<mtt,Memory_pools_use::enabled>::release_memory(void *raw_mem_ptr, std::size_t const sz) {
    if (raw_mem_ptr != nullptr) {
#ifdef MAP_MEMORY_WITH_VECTOR
        auto pos = std::find_if(data.memory_map.begin(), data.memory_map.end(),
                                [sz](std::pair <std::size_t, Memory_pool> const &x) {
                                    return x.first == sz;
                                });
#else
        auto pos = data.memory_map.find(sz);
#endif
#ifdef ALL_CHECKS_ENABLED
        if (pos == data.memory_map.end()) {
            std::cerr << "Error : Memory_mapper::release_memory(p,s) with p = @" << raw_mem_ptr
                      << " and s=" << sz << ".\n Attempting to free unallocated memory chunk.\n";
            abort();
        }
#endif
        pos->second->free(raw_mem_ptr);
    }
#ifdef ALL_CHECKS_ENABLED
    else if(sz != 0)
    {
        std::cerr << "Error : Memory_mapper::release_memory(p,s) with p = @" << raw_mem_ptr
                  << " and s=" << sz << ".\n nullptr of non null size... it seems that something has been lost (beware"
                                        " of any possible memory leak).\n";
    }
#endif
}

#endif //ifdef HAVE_BOOST

