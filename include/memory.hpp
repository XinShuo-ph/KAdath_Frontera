//
// Created by sauliac on 14/04/2020, based on a previous work from Ludwig Jens Papenfort.
//

#include "config.h"

#ifndef __MEMORY_HPP_
#define __MEMORY_HPP_


#include <array>
#include <memory>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#ifdef ENABLE_BOOST_MEMORY_POOLS
#include "boost/pool/pool.hpp"
#endif

#if MEMORY_MAP_TYPE == 3
#include "flat_hash_map.hpp"
#endif

#if MEMORY_MAP_TYPE == 0
#define MAP_MEMORY_WITH_VECTOR
#else
#define MAP_MEMORY_WITH_MAP
#endif

#define chosen_memory_map_type static_cast<Kadath::Memory::Map_type_tag>(MEMORY_MAP_TYPE)

#ifdef ENABLE_BOOST_MEMORY_POOLS
#define boost_memory_pools_use Kadath::Memory::Memory_pools_use::enabled
#else //ifdef ENABLE_BOOST_MEMORY_POOLS
#define boost_memory_pools_use Kadath::Memory::Memory_pools_use::disabled
#endif //ifdef ENABLE_BOOST_MEMORY_POOLS

namespace Kadath {
    namespace Memory {

        enum Memory_pools_use : bool {disabled = false, enabled = true};
        //! Lists the possible type to map size to already allocated memory chunk of taht size.
        enum Map_type_tag : unsigned short {
            //! value to use an implementation based on the standard library \c vector class.
            stl_vector,
            //! value to use an implementation based on the standard library \c map class.
            stl_map,
            //! value to use an implementation based on the standard library \c unoradered_map class.
            stl_unordered_map,
            //! value to use an implementation based on the Malte Skarupke hash map (see the \c flat_hash_map.hpp file).
            msc_flat_hash_map
        };

        /**
         * Struct encapsulating the memory map and pointer to memory zones list.
         * @tparam Map_type the desired map type to be used
         * @tparam memory_pools_use flag to indicate if memory pool are used instead of standard alocation/dealocation.
         */
        template<typename Map_type,Memory_pools_use memory_pools_use = Memory_pools_use::disabled>
        struct Memory_mapper_data_base;
        template<typename Map_type> struct Memory_mapper_data_base<Map_type,Memory_pools_use::disabled> {
            static_assert(std::is_default_constructible<Map_type>::value,"The map type must be default constructible.");
            //! Delete functor that will actually delete the memory chunk at the end of the execution.
            struct Standard_deleter {void operator()(void *raw_mem_ptr) {std::free(raw_mem_ptr);} };
            using raw_ptr_t = std::unique_ptr<void, Standard_deleter>;
            using ptr_list_t = std::vector<raw_ptr_t>;
            using mem_map_t = Map_type;

            //! Vector of pointers to memory zones really allocated (used or not).
            ptr_list_t ptr_list;
            //! Container mapping all available sizes with the corresponding list of already allocated (but not used)
            //! memory zones.
            mem_map_t memory_map;

            //! Default constructor.
            Memory_mapper_data_base() = default;
        };

        /**
         * Specialization for the case were memory pools are enabled.
         * @tparam Map_type underlying map type to map sizes to their corresponding memory pool.
         */
        template<typename Map_type> struct Memory_mapper_data_base<Map_type,Memory_pools_use::enabled> {
            static_assert(std::is_default_constructible<Map_type>::value,"The map type must be default constructible.");
            using mem_map_t = Map_type;
            //! Container mapping all available sizes with the corresponding memory pool.
            mem_map_t memory_map;
            //! Default constructor.
            Memory_mapper_data_base() = default;
        };

        //! Dummy template struct to map the preprocessor macro from the configure.h file to the actual map type.
        template<Map_type_tag map_type,Memory_pools_use memory_pools_use = Memory_pools_use::disabled>
        struct Memory_mapper_data {};
        template<> struct Memory_mapper_data<stl_vector,Memory_pools_use::disabled> :
                Memory_mapper_data_base<std::vector<std::pair<std::size_t, std::vector<void*>>>> {};
        template<> struct Memory_mapper_data<stl_map,Memory_pools_use::disabled> :
                Memory_mapper_data_base<std::map<std::size_t, std::vector<void*>>> {};
        template<> struct Memory_mapper_data<stl_unordered_map,Memory_pools_use::disabled> :
                Memory_mapper_data_base<std::unordered_map<std::size_t, std::vector<void*>>> {};
#if MEMORY_MAP_TYPE == 3
        // this one we implement only if explicitely requested because the corresponding map implementation spawn a lot
        // of warnings with the intel c++ compiler.
        template<> struct Memory_mapper_data<stl_unordered_map,Memory_pools_use::disabled> :
                Memory_mapper_data_base<ska::flat_hash_map<std::size_t, std::vector<void*>>> {};
#endif

#ifdef ENABLE_BOOST_MEMORY_POOLS
        using Memory_pool = std::unique_ptr<boost::pool<>>;
        template<> struct Memory_mapper_data<stl_vector,Memory_pools_use::enabled> :
                Memory_mapper_data_base<std::vector<std::pair<std::size_t,Memory_pool>>> {};
        template<> struct Memory_mapper_data<stl_map,Memory_pools_use::enabled> :
                Memory_mapper_data_base<std::map<std::size_t, Memory_pool>> {};
        template<> struct Memory_mapper_data<stl_unordered_map,Memory_pools_use::enabled> :
                Memory_mapper_data_base<std::unordered_map<std::size_t, Memory_pool>> {};
#if MEMORY_MAP_TYPE == 3
        // this one we implement only if explicitely requested because the corresponding map implementation spawn a lot
        // of warnings with the intel c++ compiler.
        template<> struct Memory_mapper_data<stl_unordered_map,Memory_pools_use::enabled> :
                Memory_mapper_data_base<ska::flat_hash_map<std::size_t, sMemory_pool>> {};
#endif
#endif //ifdef ENABLE_BOOST_MEMORY_POOLS

        /**
         * Class encapsulating a mapping beetween size and memory chunks of that size or memory pools of the size.
         * THIS IS NOT THREAD SAFE, but in our case doesn't really matter since multithreading is handled by the MPI
         * version where this memory mapper works just as well as in the sequential case.
         * @tparam map_type_tag type tag to specify the map implementation (mainly to be instantiated with the
         * preprocessor macro \c chosen_memory_map_type).
         * @tparam memory_pool_use tag to enable or diasble the use of memory pools.
         */
        template<Map_type_tag map_type_tag,Memory_pools_use memory_pool_use> class Memory_mapper;

        //! Specialization when memory pools are disabled.
        template<Map_type_tag map_type_tag> class Memory_mapper<map_type_tag,Memory_pools_use::disabled> {
        public:
            //! Alias toward the struct encapsulating internal data types.
            using Data = Memory_mapper_data<map_type_tag,Memory_pools_use::disabled>;
            //! A vector of pointer referencing memory chunk.
            using ptr_vec_t = std::vector<void *>;

        private:
            //! Data member holding the list of trully allocated memory chunks, and an object which maps a size to the
            //! unused (but allocated) memory chunks.
            static Data data;

        public:

            /**
             * Function used to request a free memory zone. Search the map in \data for the entry associated to \c sz.
             * If the mapped list of pointer is not empty, its first element is poped and returned, otherwize a new one
             * with size \c sz is allocated and referenced in \c data.ptr_list, then returned.
             * @param sz size of the desired memory zone.
             * @return a raw pointer to the memory zone of type \c void*.
             */
            static void *get_memory(std::size_t const sz);

            /**
             * Request a contiguous free memory zone for the a given number of object of a given type, calling the
             * pointer-to-\c void -returning overload.
             * @tparam T type of the objects for which the allocation is requested.
             * @param sz number of objects of type \c T.
             * @return a pointer to the memory zone.
             */
            template<typename T> static T *get_memory(std::size_t const sz) {
                return static_cast<T *>(get_memory(sz * sizeof(T)));
            }

            /**
             * Release a pointer by referenceing it into the map of free memory zones (without actual deallocation).
             * @param raw_mem_ptr the pointer to the memory zone to mark as being unused.
             * @param sz size of the memory zone.
             */
            static void release_memory(void *raw_mem_ptr, std::size_t const sz);

            /**
             * Typed version of the memory releasing method (just call the non-typed overload with the appropriate
             * size).
             * @tparam T type of the objects contained in the memory zone pointed to.
             * @param raw_mem_ptr pointer to the memory zone.
             * @param sz number of type \c T objects contained in the memory zone.
             */
            template<typename T> static void release_memory(void *raw_mem_ptr, std::size_t const sz) {
                release_memory(raw_mem_ptr, sz * sizeof(T));
            }
        };

        //! Specialization when memory pools are enabled.
        template<Map_type_tag map_type_tag> class Memory_mapper<map_type_tag,Memory_pools_use::enabled> {
        public:
            //! Alias toward the struct encapsulating internal data type.
            using Data = Memory_mapper_data<map_type_tag,Memory_pools_use::enabled>;

        private:
            //! Data member holding the memory pools mapped to their size.
            static Data data;

        public:

            /**
             * Request the memory pool with the corresponding size for a memory chunk.
             * @param sz size of the desired memory zone.
             * @return a pointer to the memory zone of type \c void*.
             */
            static void *get_memory(std::size_t const sz);

            /**
             * Request the memory pool with the corresponding size for a contiguous free memory zone for the a given
             * number of object of a given type, calling the pointer-to-\c void -returning overload.
             * @tparam T type of the objects for which the allocation is requested.
             * @param sz number of objects of type \c T.
             * @return a pointer to the memory zone.
             */
            template<typename T> static T *get_memory(std::size_t const sz) {
                return static_cast<T *>(get_memory(sz * sizeof(T)));
            }

            /**
             * Ask the proper memory pool to free the given pointer.
             * @param raw_mem_ptr the pointer to free.
             * @param sz size of the memory zone.
             */
            static void release_memory(void *raw_mem_ptr, std::size_t const sz);

            //! Typed version of the memory releasing method.
            template<typename T> static void release_memory(void *raw_mem_ptr, std::size_t const sz) {
                release_memory(raw_mem_ptr, sz * sizeof(T));
            }
        };

#include "memory_impl.hpp"

        /**
         * Struct redefining operator new and delete so that each class derived from it used the memory allocation/
         * deallocation process implementated in \c Memory_mapper (instantiated using the \c chosen_memory_map_type
         * value defined by the user when the build has been configured).
         */
        struct Memory_mapped {
            void *operator new(std::size_t sz) {return Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::get_memory(sz);}
            void operator delete(void *mem_ptr, std::size_t const sz) {
                Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::release_memory(mem_ptr, sz);
            }
            void *operator new[](std::size_t sz) {return Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::get_memory(sz);}
            void operator delete[](void *mem_ptr, std::size_t const sz) {
                Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::release_memory(mem_ptr, sz);
            }
        };

        //! Dummy struct to be used as function argument (especially constructors) to tell it we want to initialize
        //! memory.
        struct Initialize {};
        //! Tag to use to initialize memory when a function has an overload with the \c Initialize type in its
        //! signature.
        static constexpr Initialize initialize{};
        //! Dummy struct to be used as function argument (especially constructors) to tell it if we DON'T want to
        //! initialize memory.
        struct Do_not_initialize {};
        //! Tag to use to avoid memory initialization in functions with the \c Do_not_initialize type in its signature.
        static constexpr Do_not_initialize do_not_initialize{};
        //! Kind of traits to map a type to its default constructed value. By default, not usable, but a specialization
        //! can be used to map a default value to type which are not default-constructible (in term of contructors).
        template<typename T,bool has_defautl_ctor = std::is_default_constructible<T>::value>
        struct Default_initializer{};
        //! Specialization of \c Default_initializer for default-constructible types.
        template<typename T> struct Default_initializer<T,true> {static constexpr T value{};};
        //! Specialization of \c Default_initializer for pointer types.
        template<typename T> struct Default_initializer<T *,true> {static constexpr T *value{nullptr};};

        /**
         * Functor type that delete pointers with a nullity check before the actual deletion.
         * @tparam T pointer underlying type.
         */
        template<typename T> struct Safe_deleter {
            //! Alias toward the pointer type.
            using pointer = T *;
            //! Alias toward pointer to const type.
            using const_pointer = T *const;
            //! Deletion method (static version).
            static inline void apply(pointer &p) {if (p) {delete p; p = nullptr;}}
            //! Non static deletion function.
            inline void operator()(pointer &p) const { apply(p); }
            //! Deletion of a pointer to const.
            static inline void apply(const_pointer &p) {if (p) {delete p;p = nullptr;}}
            //! Non static deletion of a pointer to const.
            inline void operator()(const_pointer &p) const { apply(p); }
        };

        /**
         * Overload of the deletion functor for pointers referencing memory zones allocated using the \c new[] operator.
         * @tparam T Pointer underlying type.
         */
        template<typename T> struct Safe_deleter<T[]> {
            //! Alias toward the pointer type.
            using pointer = T *;
            //! Alias toward pointer to const type.
            using const_pointer = T *const;
            //! Deletion method (static version).
            static inline void apply(pointer &p) {if (p) {delete[] p; p = nullptr;} }
            //! Non static deletion function.
            inline void operator()(pointer &p) const { apply(p); }
            //! Deletion of a pointer to const.
            static inline void apply(const_pointer &p) {if (p) {delete[] p; /*p = nullptr;*/}}
            //! Non static deletion of a pointer to const.
            inline void operator()(const_pointer &p) const { apply(p); }
        };

        //! Delete function with nullity check.
        template<typename T> inline void safe_delete(T *&p) { Safe_deleter<T>::apply(p); }
        //! Delete function with nullity check for pointer to const.
        template<typename T> inline void safe_delete(T *const &p) { Safe_deleter<T>::apply(p); }

        /**
         * Smart pointer type for multiple allocation (i.e. arrays).
         * @tparam T data type
         * @tparam S type use for indexation (Kadath uses the \c int type).
         */
        template<typename T, typename S=int>
        class Memory_mapped_array : public Memory_mapped {
        public:
            //! Alias toward the data type.
            using value_type = T;
            //! An alias toward the supposed optimal type to be used to pass to function or return values in read-only
            //! mode.
            using value_read_only_type = typename std::conditional<(sizeof(T) < sizeof(T *)), T, T const &>::type;
            //! Alias toward the index/size type.
            using size_type = S;

        protected:
            //! Number of elements allocated.
            size_type size;
            //! Pointer to the memory zone.
            T *data;

        protected:
            /// Safety check used only from at().
            void range_check(size_type const i) const;
            //! Request for a memory zone of the same size and copy the stored value in it.
            T *duplicate_data() const;
            //! A Kadath-fashioned setter.
            size_type &set_size() { return size; }

        public:
            //! Default constructor.
            Memory_mapped_array() : size{0}, data{nullptr} {}
            //! Constructor used to convert a \c nullptr to an actual instance of \c Memory_mapped_array.
            Memory_mapped_array(std::nullptr_t) : size{0}, data{nullptr} {}
            /**
             * Main constructor : allocate a memory for the requested number of elements without initializing values.
             * @param _size number of elements.
             * @param init placeholder for overload arguments substitution.
             */
            Memory_mapped_array(size_type const _size, Do_not_initialize const init = do_not_initialize);
            /**
             * Main constructor with data values initialization.
             * @param _size number of elements
             * @param _init placeholder for overload arguments substitution.
             */
            Memory_mapped_array(size_type const _size, Initialize const _init);
            //! Copy constructor.
            Memory_mapped_array(Memory_mapped_array const &);
            //! Move constructor.
            Memory_mapped_array(Memory_mapped_array &&) noexcept;
            //! Destructor.
            virtual ~Memory_mapped_array() { this->clear(); }
            //! Copy assignment operator.
            Memory_mapped_array &operator=(Memory_mapped_array const &);
            //!  Move assignment operator.
            Memory_mapped_array &operator=(Memory_mapped_array &&) noexcept;
            //! Swaps the contents of the two instances.
            void swap(Memory_mapped_array<T> &so) noexcept {
                std::swap(size, so.size);
                std::swap(data, so.data);
            }
            /**
             * Reallocate the data and change the number of elements - contents is disacarded.
             * @param new_size
             */
            void resize(size_type const new_size);
            //! Release the memory and sets the size to zero.
            void clear() {
                Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::release_memory<T>(data, static_cast<std::size_t>(size));
                data = nullptr;
                size = 0;
            }
            /**
             * Release memory after having called a functor on each elements (usually a deleter when dealing with
             * array of arrays).
             * @tparam F the functor type.
             * @param f the functionnal object to call on each elements
             */
            template<typename F> void clear(F f) {this->apply(f); this->clear();}
            /**
             * Apply a functor to each elements of the array.
             * @tparam F functor type.
             * @param f functor to call.
             */
            template<typename F> void apply(F f) { for (auto &x : (*this)) f(x); }
            //! Returns a read-only copy of \c data (so that the class can use range-based for loop or iterator-like syntax).
            T const *begin() const noexcept { return data; };
            //! Read-only version of \c data (so that the class can use range-based for loop or iterator-like syntax).
            T const *cbegin() const noexcept { return data; };
            //! Read/write copy of \c data (so that the class can use range-based for loop or iterator-like syntax).
            T *begin() noexcept { return data; }
            //! Returns a pointer to const pointing at the end of the memory zone hold by this object.
            T const *end() const noexcept { return data + size; }
            //! Returns a pointer at the end of the memory zone hold by this object.
            T *end() noexcept { return data + size; }
            //! Returns a pointer to const pointing at the end of the memory zone hold by this object.
            T const *cend() const noexcept { return data + size; }
            //! Returns the first element (read only mode), if \c data cannot be dereference the behaviour is undefined.
            value_read_only_type front() const noexcept { return data[0]; }
            //! Returns a reference to the first element, if \c data cannot be dereference the behaviour is undefined.
            T &front() noexcept { return data[0]; }
            //! Returns the last element (read only), if \c data cannot be dereference the behaviour is undefined.
            value_read_only_type back() const noexcept { return data[size - 1]; }
            //! Returns a reference to the last element, if \c data cannot be dereference the behaviour is undefined.
            T &back() noexcept { return data[size - 1]; }
            //! Returns the value of the \c size data member in read only mode/
            size_type get_size() const noexcept { return size; }
            //! Return a pointer to const to the memory zone pointer to by \c data.
            T const *get_data() const noexcept { return data; }
            //! Kadath-fashioned setter for the \c data member.
            T *&set_data() noexcept { return data; }
            /**
             * Read-only access to elements.
             * @param i index of the elements to access to.
             * @return the value of the requested element. If \c data cannot be dereferenced or \c i is greater than or
             * equal to the value of \c size, the result is undefined.
             */
            value_read_only_type operator[](size_type const i) const noexcept { return data[i]; }
            /**
             * Read/write access to elements.
             * @param i index of the elements to access to.
             * @return a reference to the requested element. If \c data cannot be dereferenced or \c i is greater than
             * or equal to the value of \c size, the result is undefined.
             */
            T &operator[](size_type const i) noexcept { return data[i]; }
            /**
             * Read-only access to elements with safety checks.
             * @param i index of the elements to access to.
             * @return the value of the requested element. If \c data cannot be dereferenced or \c i is greater than or
             * equal to the value of \c size, a runtime error exception is thrown.
             */
            value_read_only_type at(size_type const i) const {this->range_check(i);return data[i];}
            /**
             * Read/write access to elements.
             * @param i index of the elements to access to.
             * @return a reference to the requested element. If \c data cannot be dereferenced or \c i is greater than or
             * equal to the value of \c size, a runtime error exception is thrown.
             */
            T &at(size_type const i) {this->range_check(i);return data[i];}
            //! Emptyness check.
            bool empty() const noexcept { return size == 0; }
        };

        template<typename T, typename S>
        inline Memory_mapped_array<T, S>::Memory_mapped_array(size_type const _size, Do_not_initialize const)
            : size{_size}, data{Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::get_memory<T>(static_cast<std::size_t>(size))} {}

        template<typename T, typename S>
        inline Memory_mapped_array<T, S>::Memory_mapped_array(size_type const _size, Initialize const)
            : size{_size}, data{Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::get_memory<T>(static_cast<std::size_t>(size))}
        {
            for (auto &x:(*this)) x = Default_initializer<T>::value;
        }

        template<typename T, typename S>
        inline Memory_mapped_array<T, S>::Memory_mapped_array(Memory_mapped_array const &source)
                : size{source.size}, data{source.duplicate_data()} {}

        template<typename T, typename S>
        inline Memory_mapped_array<T, S>::Memory_mapped_array(Memory_mapped_array &&source) noexcept
                : size{source.size}, data{source.data} {
            source.size = 0;
            source.data = nullptr;
        }

        template<typename T, typename S>
        inline Memory_mapped_array<T, S> &
        Memory_mapped_array<T, S>::operator=(Memory_mapped_array<T, S> const &source) {
            this->resize(source.size);
            for (std::size_t i{0}; i < size; i++) data[i] = source.data[i];
            return *this;
        }

        template<typename T, typename S>
        inline Memory_mapped_array<T, S> &
        Memory_mapped_array<T, S>::operator=(Memory_mapped_array<T, S> &&source) noexcept { this->swap(source); }

        template<typename T, typename S>
        void Memory_mapped_array<T, S>::range_check(const size_type i) const {
            if (i >= size || data == nullptr) {
                throw std::runtime_error{
                        std::string{"Memory_mapped_array::range_check : i (which is " + std::to_string(i) +
                                    ") >= this->set_size() (which is " + std::to_string(size) + ")."}};
            }
        }

        template<typename T, typename S>
        void Memory_mapped_array<T, S>::resize(size_type const new_size) {
            if (new_size != size) {
                Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::release_memory<T>(data, size);
                size = new_size;
                data = Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::get_memory<T>(static_cast<std::size_t>(size));
            }
        }

        template<typename T, typename S>
        inline T *Memory_mapped_array<T, S>::duplicate_data() const {
            T *const data_copy{Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::get_memory<T>(static_cast<std::size_t>(size))};
            for (size_type i{0}; i < size; i++) data_copy[i] = data[i];
            return data_copy;
        }

        /**
         * An allocator type in the STL fashion, so that STL container may use the \c Memory_mapper class.
         * @tparam T type of the elements.
         * @tparam S size/index type.
         */
        template<typename T, typename S=std::size_t> struct Memory_mapped_allocator {
            using value_type = T;
            using size_type = S;
            using pointer = T *;
            using const_pointer = T const *;
            using reference = T &;
            using const_reference = T const &;
            using difference_type = std::ptrdiff_t;
            using propagate_on_container_move_assignment = std::true_type;
            using is_always_equal = std::true_type;

            //! Default constructor.
            Memory_mapped_allocator() = default;
            //! Copy constructor.
            template<typename U> constexpr Memory_mapped_allocator(const Memory_mapped_allocator<U> &) noexcept {}
            pointer address(reference x) const noexcept { return &x; }
            const_pointer address(const_reference x) const noexcept { return &x; }
            pointer allocate(size_type sz) { return Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::get_memory<T>(static_cast<std::size_t>(sz)); }
            void deallocate(T *p, size_type sz) noexcept {
                Memory_mapper<chosen_memory_map_type,boost_memory_pools_use>::release_memory<T>(p, static_cast<std::size_t>(sz));
            }

        };
        //! Checks wether allocator for different element types can use the same memory mapper (always true).
        template<class T, class U>
        bool
        operator==(const Memory_mapped_allocator<T> &, const Memory_mapped_allocator<U> &) noexcept { return true; }
        //! Checks if allocator for different element types cannot use the same memory mapper (always false).
        template<class T, class U>
        bool
        operator!=(const Memory_mapped_allocator<T> &, const Memory_mapped_allocator<U> &) noexcept { return false; }

        //! An alias that somewhat change the default allocator of the STL vector and use our \c Memory_mapped_allocator.
        template<typename T, typename S> using Memory_mapped_vector = std::vector<T, Memory_mapped_allocator<T, S>>;



        template<typename T> using MMPtr_array = Memory_mapped_array<T *>;
    }
    // To get back the needed symbols in the \c Kadath namespace.
    using Memory::Memory_mapper;
    using Memory::Memory_mapped;
    using Memory::Memory_mapped_array;
    using Memory::Memory_mapped_allocator;
    using Memory::Safe_deleter;
    using Memory::safe_delete;
    using Memory::Initialize;
    using Memory::initialize;
    using Memory::Do_not_initialize;
    using Memory::do_not_initialize;
    using Memory::MMPtr_array;
    using Memory::Memory_mapped_vector;
}


#endif //__MEMORY_HPP_