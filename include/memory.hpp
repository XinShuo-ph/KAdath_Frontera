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
#else
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
        size_type size;
        T* data;

    protected:
        /// Safety check used only from at().
        void range_check(size_type const i) const;
        T* duplicate_data() const;
        size_type & set_size() {return size;}

    public:
        Memory_mapped_array() : size{0}, data{nullptr} {}
        Memory_mapped_array(std::nullptr_t ) : size{0}, data{nullptr} {}
        Memory_mapped_array(Initialize const) : size{0}, data{nullptr} {for(auto & x : (*this)) x = T{};}
        Memory_mapped_array(size_type const,Do_not_initialize const init = do_not_initialize);
        Memory_mapped_array(size_type const,Initialize const);
        Memory_mapped_array(Memory_mapped_array const &);
        Memory_mapped_array(Memory_mapped_array &&);
        virtual ~Memory_mapped_array() {this->clear();}

        Memory_mapped_array & operator=(Memory_mapped_array const &);
        Memory_mapped_array & operator=(Memory_mapped_array &&);
        void swap(Memory_mapped_array<T> &so) noexcept {std::swap(size,so.size); std::swap(data,so.data);}
        void resize(size_type const new_size);
        void clear() {Memory_mapper::release_memory<T>(data,static_cast<std::size_t>(size)); data=nullptr; size=0;}
        template<typename F> void clear(F f) {this->apply(f); this->clear();}
        template<typename F> void apply(F f) {for(auto &x : (*this)) f(x);}

        T const * begin() const noexcept {return data;};
        T const * cbegin() const noexcept {return data;};
        T * begin() noexcept {return data;}
        T const * end() const noexcept {return data + size;}
        T * end() noexcept {return data + size;}
        T const * cend() const noexcept {return data + size;}
        value_read_only_type front() const noexcept {return data[0];}
        T & front() noexcept {return data[0];}
        value_read_only_type back() const noexcept {return data[size-1];}
        T & back() noexcept {return data[size-1];}

        size_type get_size() const noexcept {return size;}
        T const * get_data() const noexcept {return data;}
        T*& set_data() noexcept {return data;}
        value_read_only_type operator[](size_type const i) const noexcept {return data[i];}
        T& operator[](size_type const i) noexcept {return data[i];}
        value_read_only_type at(size_type const i) const {this->range_check(i); return data[i];}
        T& at(size_type const i) {this->range_check(i); return data[i];}

        bool empty() const noexcept {return size==0;}

//        operator T const * () const {return data;}
//        operator T * () {return data;}
    };

    template<typename T,typename S> inline Memory_mapped_array<T,S>::Memory_mapped_array(size_type const _size,
                                                                            Do_not_initialize const)
            : size{_size}, data{Memory_mapper::get_memory<T>(static_cast<std::size_t>(size))} {}

    template<typename T,typename S> inline Memory_mapped_array<T,S>::Memory_mapped_array(size_type const _size,Initialize const)
            : size{_size}, data{Memory_mapper::get_memory<T>(static_cast<std::size_t>(size))} {for(auto & x:(*this)) x = T{};}

    template<typename T,typename S> inline Memory_mapped_array<T,S>::Memory_mapped_array(Memory_mapped_array const & source)
            : size{source.size}, data{source.duplicate_data()} {}
    template<typename T,typename S> inline Memory_mapped_array<T,S>::Memory_mapped_array(Memory_mapped_array && source)
            : size{source.size}, data{source.data} {source.size = 0; source.data = nullptr;}

    template<typename T,typename S> inline Memory_mapped_array<T,S> &
            Memory_mapped_array<T,S>::operator=(Memory_mapped_array<T,S> const &source) {
                if(size != source.size) {
                    Memory_mapper::release_memory<T>(data,static_cast<std::size_t>(size));
                    size = source.size;
                }
                for(std::size_t i{0};i<size;i++) data[i] = source.data[i];
                return *this;
            }
    template<typename T,typename S> inline Memory_mapped_array<T,S> &
            Memory_mapped_array<T,S>::operator=(Memory_mapped_array<T,S> && source) {this->swap(source);}

    template<typename T,typename S> void Memory_mapped_array<T,S>::range_check(const size_type i) const {
        if (i >= this->set_size()) {
            throw std::runtime_error{
                    std::string{"Memory_mapped_array::range_check : i (which is " + std::to_string(i) +
                                ") >= this->set_size() (which is " + std::to_string(size) + ")."}};
        }
    }

    template<typename T,typename S> void Memory_mapped_array<T,S>::resize(size_type const new_size) {
        if(new_size != size) {
            Memory_mapper::release_memory<T>(data,size);
            size = new_size;
            data = Memory_mapper::get_memory<T>(static_cast<std::size_t>(size));
        }
    }

    template<typename T,typename S> inline T * Memory_mapped_array<T,S>::duplicate_data() const {
        T * const data_copy {Memory_mapper::get_memory<T>(static_cast<std::size_t>(size))};
        for(size_type i{0};i<size;i++) data_copy[i] = data[i];
        return data_copy;
    }


    template <typename T,typename S=std::size_t> struct Memory_mapped_allocator
    {
        using value_type = T;
        using size_type = S;
        using pointer = T*;
        using const_pointer = T const *;
        using reference = T&;
        using const_reference = T const &;
        using difference_type = std::ptrdiff_t ;
        using propagate_on_container_move_assignment = std::true_type;
        using is_always_equal = std::true_type;


        Memory_mapped_allocator () = default;
        template <typename U> constexpr Memory_mapped_allocator (const Memory_mapped_allocator <U>&) noexcept {}

        pointer address(reference x) const noexcept {return &x;}
        const_pointer address(const_reference x) const noexcept {return &x;}

        pointer allocate(size_type sz) {return Memory_mapper::get_memory<T>(static_cast<std::size_t>(sz));}
        void deallocate(T* p, size_type sz) noexcept { Memory_mapper::release_memory<T>(p,static_cast<std::size_t>(sz)); }

    };

    template <class T, class U>
    bool operator==(const Memory_mapped_allocator <T>&, const Memory_mapped_allocator <U>&) noexcept { return true; }
    template <class T, class U>
    bool operator!=(const Memory_mapped_allocator <T>&, const Memory_mapped_allocator <U>&) noexcept { return false; }

    template<typename T,typename S> using Memory_mapped_vector = std::vector<T,Memory_mapped_allocator<T,S>>;
}


#endif //__MEMORY_HPP_