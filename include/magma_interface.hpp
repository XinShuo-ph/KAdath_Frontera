/*
    Copyright 2020 sauliac

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
#ifndef __MAGMA_INTERFACE_HPP_
#define __MAGMA_INTERFACE_HPP_

#include <vector>
#include <memory>
#include "config.h"
#include "array.hpp"

#ifdef ENABLE_GPU_USE
#include "magma_v2.h"
#include "magma_lapack.h"
#endif


#define TESTING_CHECK( err )                                                 \
    do {                                                                     \
        magma_int_t err_ = (err);                                            \
        if ( err_ != 0 ) {                                                   \
            fprintf( stderr, "Error: %s\nfailed at %s:%d: error %lld: %s\n", \
                     #err, __FILE__, __LINE__,                               \
                     (long long) err_, magma_strerror(err_) );               \
            exit(1);                                                         \
        }                                                                    \
    } while( 0 )


namespace Kadath {


#ifdef ENABLE_GPU_USE
    //! An allocator to be passed as template parameter to the STL vector to assure optimal data alignment for Magma.
    template<typename T> struct Magma_allocator //: public std::allocator_traits<std::allocator<double>>
    {
        //! Used to build array of doubles.
        using value_type = T;
        using pointer = value_type *;
        using const_pointer = value_type const *;
        using void_pointer = void*;
        using const_void_pointer = void const *;
        using size_type = std::size_t;
        //! Default constructor.
        Magma_allocator() = default;
        template<class U> constexpr Magma_allocator(Magma_allocator<U> const &) noexcept {};

        pointer allocate(std::size_t n);

        void deallocate(pointer data,std::size_t) noexcept {
            magma_free_cpu(data);
        }
    };

    template<> double * Magma_allocator<double>::allocate(std::size_t n)
    {
        double * data{nullptr};
        magma_int_t err{magma_dmalloc_cpu(&data,static_cast<magma_int_t>(n))};
        if(err!= 0) throw std::bad_alloc{};
        return data;
    }
    template<> magma_int_t * Magma_allocator<magma_int_t>::allocate(std::size_t n)
    {
        magma_int_t * data{nullptr};
        magma_int_t err{magma_imalloc_cpu(&data,static_cast<magma_int_t>(n))};
        if(err!= 0) throw std::bad_alloc{};
        return data;
    }

    //for the moment, no rebind only magma double allocators
    template<typename A,typename B> inline constexpr bool operator==(Magma_allocator<A> const &,Magma_allocator<B> const &)
    {return std::is_same<A,B>::value;}
    template<typename A,typename B> inline constexpr bool operator!=(Magma_allocator<A> const &,Magma_allocator<B> const &)
    {return !std::is_same<A,B>::value;}


    class Magma_array : public std::vector<double,Magma_allocator<double>>
    {
    protected:
        using Base = std::vector<double,Magma_allocator<double>>;
        magma_int_t dim;
    public:
        Magma_array(std::size_t _dim) :  Base(_dim),dim{static_cast<magma_int_t>(_dim)} {}
        Magma_array(Array<double> const & source) ;
        magma_int_t get_dim() const {return dim;}
        //! Copy the source array to \c this without reallocation so the sizes must match.
        Magma_array & operator=(Array<double> const &source);
    };
    /**
     * Class encapsulating an array of double representing a square matrix to be manipulated with magma linear solving
     * routines.
     */
    class Magma_matrix : public Magma_array
    {
    public:
        using PPivot = std::unique_ptr<std::vector<magma_int_t,Magma_allocator<magma_int_t>>>;
    private:
        //! Order of the matrix.
        magma_int_t order;
        //! Leading dimension array.
        magma_int_t lda;
        //!
        PPivot pivot;

    public:
        //! Allocation of data without initialization.
        Magma_matrix(std::size_t _order,std::size_t _lda=0) :
            Magma_array{_lda==0 ? _order*_order : _order*_lda},
            order{static_cast<magma_int_t>(_order)},
            lda{static_cast<magma_int_t>(_lda==0 ? _order : _lda)}, pivot{nullptr}
        {assert(lda >= std::max(1,order));}
        //! Allocate data and initializes it using the values of the source array (must be a column major index matrix).
        Magma_matrix(Array<double> const & source,std::size_t n) : Magma_array{source},
                                                                   order{static_cast<magma_int_t>(n)},
                                                                   lda{static_cast<magma_int_t>(n)}, pivot{nullptr} {}

        bool is_factorized() const {return pivot != nullptr;}
        double operator()(std::size_t i,std::size_t j) const
        {
            assert(0<=i && i<dim && 0<=j && j<lda);
            return (*this)[i+lda*j];
        }
        //! Writable row/col access to data.
        double & operator()(std::size_t i,std::size_t j)
        {
            assert(0<=i && i<dim && 0<=j && j<lda);
            return (*this)[i+lda*j];
        }

        Magma_matrix & set_column(std::size_t col,Array<double> const & source)
        {
            assert(source.get_nbr() == order);
            for(std::size_t i{0};i<source.get_nbr();i++)
            {
                (*this)(i,col) = source.get_data()[i];
            }
            return *this;
        }

        magma_int_t get_order() const {return order;}
        magma_int_t get_lda() const {return lda;}

        Magma_array & solve(Magma_array & second_member);
	using Magma_array::operator=;
    };

#endif

}
#endif //__MAGMA_INTERFACE_HPP_
