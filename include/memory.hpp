//
// Created by sauliac on 14/04/2020.
//

#include <array>
#include <memory>
#include <vector>
#include <algorithm>

#include "config.h"

#if MEMORY_MAP_TYPE == 1
#include <unordered_map>
#elif MEMORY_MAP_TYPE == 2
#include <map>
#elif MEMORY_MAP_TYPE == 3
#include "flat_hash_map.hpp"
#endif

#if MEMORY_MAP_TYPE == 0
#define MAP_MEMORY_WITH_VECTOR
#else
#define MAP_MEMORY_WITH_MAP
#endif


#ifndef __MEMORY_HPP_
#define __MEMORY_HPP_

namespace Kadath {

    class Coef_mem {
    private:
        static std::array<std::unique_ptr<double[]>,2> mem_ptrs;
        static std::array<std::size_t,2> lengths;

    public:
        static double* get_mem(std::size_t const num, std::size_t const len) {
            if(len > lengths[num]) {
                mem_ptrs[num].reset(new double[len]);
                lengths[num] = len;
            }
            return mem_ptrs[num].get();
        }
    };


    class Memory_mapper {
    private:
        struct Free_deleter {
            void operator()(void* raw_mem_ptr) { std::free(raw_mem_ptr); }
        };

        using ptr_vec_t = std::vector<void*>;

#if MEMORY_MAP_TYPE == 0
        using mem_map_t = std::vector<ptr_vec_t>;
        using mem_sizes_t = std::vector<std::size_t>;
        static mem_sizes_t mem_sizes;
#elif MEMORY_MAP_TYPE == 1
        using mem_map_t = std::unordered_map<std::size_t,ptr_vec_t>;
#elif MEMORY_MAP_TYPE == 2
        using mem_map_t = std::map<std::size_t,ptr_vec_t>;
#elif MEMORY_MAP_TYPE == 3
        using mem_map_t = ska::flat_hash_map<std::size_t, ptr_vec_t>;

#endif
        using raw_ptr_t = std::unique_ptr<void,Free_deleter>;
        using ptr_list_t = std::vector<raw_ptr_t>;


        static mem_map_t memory_map;
        static ptr_list_t ptr_list;

    public:

        static void* get_memory(std::size_t const sz) {
            if(sz == 0) return nullptr;
            void* raw_mem_ptr;
#ifdef MAP_MEMORY_WITH_VECTOR
            auto pos = std::find(std::begin(mem_sizes), std::end(mem_sizes), sz);
            if(pos == mem_sizes.end()) {
                mem_sizes.push_back(sz);
                memory_map.resize(mem_sizes.size());
                pos = --mem_sizes.end();
            }

            std::size_t index = pos - mem_sizes.begin();
            auto& mem = memory_map[index];
#else
            auto& mem = memory_map[sz];
#endif
            if(mem.empty()) {
                ptr_list.emplace_back(std::malloc(sz));
                raw_mem_ptr = ptr_list.back().get();
            }
            else {
                raw_mem_ptr = mem.back();
                mem.pop_back();
            }

            return raw_mem_ptr;
        }

        template<typename T>
        static T* get_memory(std::size_t const sz) {
            return static_cast<T*>(get_memory(sz * sizeof(T)));
        }

        static void release_memory(void* raw_mem_ptr, std::size_t const sz)  {
            if(raw_mem_ptr != nullptr) {
#if MEMORY_MAP_TYPE == 0
                auto pos = std::find(std::begin(mem_sizes), std::end(mem_sizes), sz);
                std::size_t index = pos - mem_sizes.begin();
                memory_map[index].push_back(raw_mem_ptr);
#elif MEMORY_MAP_TYPE == 3
                memory_map[sz].push_back(raw_mem_ptr);
#endif
            }
        }

        template<typename T>
        static void release_memory(void* raw_mem_ptr, std::size_t const sz) {
            release_memory(raw_mem_ptr, sz * sizeof(T));
        }
    };

    struct Memory_mapped {
        void* operator new(std::size_t sz) {
            return Memory_mapper::get_memory(sz);
        }

        void operator delete(void* mem_ptr, std::size_t const sz) {
            Memory_mapper::release_memory(mem_ptr, sz);
        }

        void* operator new[](std::size_t sz) {
            return Memory_mapper::get_memory(sz);
        }

        void operator delete[](void* mem_ptr, std::size_t const sz) {
            Memory_mapper::release_memory(mem_ptr, sz);
        }
    };

}

void* operator new(std::size_t sz) {return Kadath::Memory_mapper::get_memory(sz);}

void operator delete(void* mem_ptr, std::size_t const sz) {Kadath::Memory_mapper::release_memory(mem_ptr, sz);}

void* operator new[](std::size_t sz) {return Kadath::Memory_mapper::get_memory(sz);}

void operator delete[](void* mem_ptr, std::size_t const sz) {Kadath::Memory_mapper::release_memory(mem_ptr, sz);}

#endif //__MEMORY_HPP_
