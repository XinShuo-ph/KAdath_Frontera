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
#include <config.h>

#ifdef PAR_VERSION
#include "mpi.h"
#endif

#include <vector>
#include <type_traits>

#ifdef ENABLE_GPU_USE
#include "magma_v2.h"
#include "magma_lapack.h"
#endif

#include "system_of_eqs.hpp"
#include "matrice.hpp"
#include "scalar.hpp"
#include "array_math.cpp"


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

        double * allocate(std::size_t n);

        void deallocate(double * data,std::size_t) noexcept;
    };

    template<> double * Magma_allocator<double>::allocate(std::size_t n)
    {
        double * data{nullptr};
        magma_int_t err{magma_dmalloc_cpu(&data,static_cast<magma_int_t>(n))};
        if(err!= 0) throw std::bad_alloc{};
        return data;
    }
    template<> void Magma_allocator<double>::deallocate(double *data, std::size_t) noexcept
    {
        magma_free_cpu(data);
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
    private:
        //! Order of the matrix.
        magma_int_t order;
        //! Leading dimension array.
        magma_int_t lda;

    public:
        //! Allocation of data without initialization.
        Magma_matrix(std::size_t _order,std::size_t _lda=0) :
            Magma_array{_lda==0 ? _order*_order : _order*_lda},
            order{static_cast<magma_int_t>(_order)},
            lda{static_cast<magma_int_t>(_lda==0 ? _order : _lda)}
        {assert(lda >= std::max(1,order));}
        //! Allocate data and initializes it using the values of the source array (must be a column major index matrix).
        Magma_matrix(Array<double> const & source,std::size_t n) : Magma_array{source},
                                                                   order{static_cast<magma_int_t>(n)},
                                                                   lda{static_cast<magma_int_t>(n)} {}

        bool is_factorized() const {return false;}
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

        magma_int_t get_order() const {return order;}
        magma_int_t get_lda() const {return lda;}

    };

#endif

    template<>
    bool System_of_eqs::do_newton<Computational_model::gpu_sequential>(double precision, double& error,std::ostream &os,
                                                                   Array<double> * copy_matrix)
    {
#ifdef PAR_VERSION
        int rank;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        if(rank==0) {
#endif
            niter++;
            if(niter==1 && display_newton_data)
            {
                display_do_newton_report_header(os,precision);
            }
            Array<double> second(sec_member());
            error = max(fabs(second));
            if (error < precision)
            {
                if(display_newton_data)
                {
                    display_do_newton_ending_line(os,precision,error);
                    os  << endl;
                }
                return true;
            } else {
                int nn(second.get_size(0));
                if (nbr_unknowns != nn) {
                    cerr << "N unknowns  = " << nbr_unknowns << endl;
                    cerr << "N equations = " << nn << endl;
                    abort();
                }

                Hash_key chrono_key = this->start_chrono("do_newton | problem size = ", nn, " | matrix computation ");
                Matrice ope(nn, nn);
                compute_matrix_adjacent(ope.get_array(), nn);
#ifdef ENABLE_GPU_USE
                Magma_matrix magma_matrix{ope.get_array(),static_cast<std::size_t>(nn)};
#endif
                if(copy_matrix && niter==1) *copy_matrix = ope.get_array();
                Duration const
                        t_load_matrix{this->stop_chrono(chrono_key)};


                chrono_key = this->start_chrono("do_newton | problem size = ", nn, " | matrix inversion ");
                ope.set_lu();
                Array<double> xx(ope.solve(second));
                Duration const
                        t_inv_matrix{this->stop_chrono(chrono_key)};

                chrono_key = this->start_chrono("do_newton | problem size = ", nn, " | Newton update ");
                newton_update_vars(xx);
                Duration const t_newton_update
                        {this->stop_chrono(chrono_key)};
                if(display_newton_data)
                {
                    display_do_newton_iteration(os,
                                                {niter,nn,error,t_load_matrix,Duration{},t_inv_matrix,t_newton_update});
                }
                return false;
            }
#ifdef PAR_VERSION
        }
        else {
            return false;
        }
#endif
    }
//#endif

#ifdef ENABLE_GPU_USE
    Magma_array::Magma_array(const Array<double> &source) : Base(source.get_nbr()),dim{static_cast<magma_int_t>(source.get_nbr())}
    {
#ifdef PAR_VERSION
        int rank, nproc;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        MPI_Comm_size (MPI_COMM_WORLD, &nproc);
        std::size_t const r {this->size() % static_cast<std::size_t>(nproc)};
        std::size_t const nb_coeff{this->size() / nproc + ((rank < r) ? 1:0)};
        std::size_t const i_start{rank * nb_coeff + (rank < r ? static_cast<std::size_t>(rank) : r)};
#else
        std::size_t const i_start{0};
        std::size_t const nb_coeff{this->size()};
#endif
        for(std::size_t i{i_start};i<nb_coeff;i++)
        {
            std::size_t const k{i+i_start};
            (*this)[k] = source.get_data()[k];
        }
    }

    Magma_array & Magma_array::operator=(const Array<double> &source)
    {
        assert(this->size() == source.get_nbr());
#ifdef PAR_VERSION
        int rank, nproc;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        MPI_Comm_size (MPI_COMM_WORLD, &nproc);
        std::size_t const r {this->size() % static_cast<std::size_t>(nproc)};
        std::size_t const nb_coeff{this->size() / nproc + ((rank < r) ? 1:0)};
        std::size_t const i_start{rank * nb_coeff + (rank < r ? static_cast<std::size_t>(rank) : r)};
#else
        std::size_t const i_start{0};
        std::size_t const nb_coeff{this->size()};
#endif
        for(std::size_t i{i_start};i<nb_coeff;i++)
        {
            std::size_t const k{i+i_start};
            (*this)[k] = source.get_data()[k];
        }
        return *this;
    }
#endif

}

