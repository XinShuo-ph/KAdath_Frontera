//
// Created by sauliac on 14/04/2020.
//

#include "config.h"

#ifndef __MEMORY_HPP_
#define __MEMORY_HPP_


#include <array>
#include <memory>
#include <vector>

#include <algorithm>

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

namespace Kadath {

    class Coef_mem {
    private:
        static std::array<std::unique_ptr<double[]>, 2> mem_ptrs;
        static std::array<std::size_t, 2> lengths;

    public:
        static double *get_mem(std::size_t const num, std::size_t const len) {
            if (len > lengths[num]) {
                mem_ptrs[num].reset(new double[len]);
                lengths[num] = len;
            }
            return mem_ptrs[num].get();
        }
    };


    class Memory_mapper {
    private:
        struct Free_deleter {
            void operator()(void *raw_mem_ptr) { std::free(raw_mem_ptr); }
        };

        using ptr_vec_t = std::vector<void *>;

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
        using raw_ptr_t = std::unique_ptr<void, Free_deleter>;
        using ptr_list_t = std::vector<raw_ptr_t>;


        static mem_map_t memory_map;
        static ptr_list_t ptr_list;

    public:

        static void *get_memory(std::size_t const sz) {
            if (sz == 0) return nullptr;
            else {
                void *raw_mem_ptr;
#ifdef MAP_MEMORY_WITH_VECTOR
                auto pos = std::find(std::begin(mem_sizes), std::end(mem_sizes), sz);
                if (pos == mem_sizes.end()) {
                    mem_sizes.push_back(sz);
                    memory_map.resize(mem_sizes.size());
                    pos = --mem_sizes.end();
                }

                std::size_t index = pos - mem_sizes.begin();
                auto &mem = memory_map[index];
#else
                auto& mem = memory_map[sz];
#endif
                if (mem.empty()) {
                    ptr_list.emplace_back(std::malloc(sz));
                    raw_mem_ptr = ptr_list.back().get();
                } else {
                    raw_mem_ptr = mem.back();
                    mem.pop_back();
                }

                return raw_mem_ptr;
            }
        }

        template<typename T>
        static T *get_memory(std::size_t const sz) {
            return static_cast<T *>(get_memory(sz * sizeof(T)));
        }

        static void release_memory(void *raw_mem_ptr, std::size_t const sz) {
            if (raw_mem_ptr != nullptr) {
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
        static void release_memory(void *raw_mem_ptr, std::size_t const sz) {
            release_memory(raw_mem_ptr, sz * sizeof(T));
        }
    };

    struct Memory_mapped {
        void *operator new(std::size_t sz) {
            return Memory_mapper::get_memory(sz);
        }

        void operator delete(void *mem_ptr, std::size_t const sz) {
            Memory_mapper::release_memory(mem_ptr, sz);
        }

        void *operator new[](std::size_t sz) {
            return Memory_mapper::get_memory(sz);
        }

        void operator delete[](void *mem_ptr, std::size_t const sz) {
            Memory_mapper::release_memory(mem_ptr, sz);
        }
    };


    template<typename T,typename S=int> class Memory_mapped_array : public Memory_mapped {
    public:
        struct Initialize {};
        static constexpr Initialize initialize{};
        struct Do_not_initialize {};
        static constexpr Do_not_initialize do_not_initialize{};
        using value_type = T;
        using value_read_only_type = typename std::conditional< (sizeof(T)<sizeof(T*)),T,T const &>::type;
        using size_type = S;

    protected:
        size_type m_size;
        T* m_data;

    protected:
        /// Safety check used only from at().
        void range_check(size_type const i) const;
        T* duplicate_data() const;
        size_type & size() {return m_size;}
        T*& data() noexcept {return m_data;}

    public:
        Memory_mapped_array() : m_size{0}, m_data{nullptr} {}
        Memory_mapped_array(Initialize const) : m_size{0}, m_data{nullptr} {for(auto & x : (*this)) x = T{};}
        Memory_mapped_array(size_type const,Do_not_initialize const init = do_not_initialize);
        Memory_mapped_array(size_type const,Initialize const);
        Memory_mapped_array(Memory_mapped_array const &);
        Memory_mapped_array(Memory_mapped_array &&);
        virtual ~Memory_mapped_array() {this->clear();}

        Memory_mapped_array & operator=(Memory_mapped_array const &) = delete;
        Memory_mapped_array & operator=(Memory_mapped_array &&);
        void swap(Memory_mapped_array<T> &so) noexcept {std::swap(m_size,so.m_size); std::swap(m_data,so.m_data);}
        void resize(size_type const new_size);
        void clear() {Memory_mapper::release_memory<T>(m_data,static_cast<std::size_t>(m_size)); m_data=nullptr; m_size=0;}

        T const * begin() const noexcept {return data;};
        T const * cbegin() const noexcept {return data;};
        T * begin() noexcept {return data;}
        T const * end() const noexcept {return m_data + m_size;}
        T * end() noexcept {return m_data + m_size;}
        T const * cend() const noexcept {return m_data + m_size;}
        value_read_only_type front() const noexcept {return m_data[0];}
        T & front() noexcept {return m_data[0];}
        value_read_only_type back() const noexcept {return m_data[m_size-1];}
        T & back() noexcept {return m_data[m_size-1];}

        size_type size() const noexcept {return m_size;}
        T const * data() const noexcept {return m_data;}
        value_read_only_type operator[](size_type const i) const noexcept {return m_data[i];}
        T& operator[](size_type const i) noexcept {return m_data[i];}
        value_read_only_type at(size_type const i) const {this->range_check(i); return m_data[i];}
        T& at(size_type const i) {this->range_check(i); return m_data[i];}
    };

    template<typename T,typename S> inline Memory_mapped_array<T,S>::Memory_mapped_array(size_type const _size,
                                                                            Do_not_initialize const)
            : m_size{_size}, m_data{Memory_mapper::get_memory<T>(static_cast<std::size_t>(size))} {}

    template<typename T,typename S> inline Memory_mapped_array<T,S>::Memory_mapped_array(size_type const _size,Initialize const)
            : m_size{_size}, m_data{Memory_mapper::get_memory<T>(static_cast<std::size_t>(size))} {for(auto & x:(*this)) x = T{};}

    template<typename T,typename S> inline Memory_mapped_array<T,S>::Memory_mapped_array(Memory_mapped_array const & source)
            : m_size{source.size}, m_data{source.duplicate_data()} {}
    template<typename T,typename S> inline Memory_mapped_array<T,S>::Memory_mapped_array(Memory_mapped_array && source)
            : m_size{source.size}, m_data{source.data} {source.size = 0; source.data = nullptr;}

    template<typename T,typename S> inline Memory_mapped_array<T,S> &
            Memory_mapped_array<T,S>::operator=(Memory_mapped_array<T,S> && source) {this->swap(source);}

    template<typename T,typename S> void Memory_mapped_array<T,S>::range_check(const size_type i) const {
        if (i >= this->size()) {
            throw std::runtime_error{
                    std::string{"Memory_mapped_array::range_check : i (which is " + std::to_string(i) +
                                ") >= this->size() (which is " + std::to_string(m_size) + ")."}};
        }
    }

    template<typename T,typename S> void Memory_mapped_array<T,S>::resize(size_type const new_size) {
        if(new_size != m_size) {
            Memory_mapper::release_memory<T>(m_data,m_size);
            m_size = new_size;
            m_data = Memory_mapper::get_memory<T>(static_cast<std::size_t>(m_size));
        }
    }

    template<typename T,typename S> inline T * Memory_mapped_array<T,S>::duplicate_data() const {
        T * const data_copy {Memory_mapper::get_memory<T>(static_cast<std::size_t>(m_size))};
        for(size_type i{0};i<m_size;i++) data_copy[i] = m_data[i];
        return data_copy;
    }

#ifdef ENABLE_CUSTOM_MEMORY_MAPPING
#define MEMORY_ALLOCATION_POLICY : Memory_mapped
#define ARRAY_ALLOCATION_POLICY(value_type) Memory_mapped_array<value_type>
#else
#define MEMORY_ALLOCATION_POLICY
#define ARRAY_ALLOCATION_POLICY(value_type) value_type *
#endif
}


#endif //__MEMORY_HPP_