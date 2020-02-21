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
#ifndef __MULTI_THREAD_HPP_
#define __MULTI_THREAD_HPP_

#include <thread>
#include "system_of_eqs.hpp"

namespace Kadath
{
    /**
     * This class is a thread-safe version of \a System_of_eqs. Thread-safety is achieved by brute-force as is done
     * with MPI by maintaining a specific copy of the original instance for each possible thread of execution. Each
     * setter method is spread toward the thread specific copies. In the case were the max number of thread is only one,
     * the class exactly behave the same as \a System_of_eqs.
     * Be aware that the class hides non-virtual inherited methods, so this should not be used through a pointer to
     * base class, unless you want only the main-thread-related instance of \a System_of_eqs to be used.
     */
    class System_of_eqs_threaded : public System_of_eqs
    {
    public:
        using container = std::vector<std::unique_ptr<System_of_eqs>>;
    private:
        std::size_t max_num_threads;
        container thread_specific_systems;

        /**
         * Function that spread a call to some method of the base class \a System_od_eqs to each instance of the
         * \a thread_specific_systems vector.
         * @tparam Args types of the arguments of the method.
         * @param method pointer to member function to be called.
         * @param args arguments values.
         */
        template<typename... Args> inline void call_base_method(void (System_of_eqs::*method)(Args...),Args... args);

    public:
        /**
         * Template constructor transfering all constructors from \a System_of_eqs to this class.
         * @tparam Args any template parameter pack with same types as any of one of \a System_of_eqs constructors.
         * @param max_num_threads maximum number of threads (use the number of cors of each nodes).
         * @param args arguments values of the template parameter pack.
         */
        template<typename... Args>
        System_of_eqs_threaded(std::size_t max_num_threads,Args && ... args);

        /**
         * Multi-threaded version of the \c build_matrix method.
         * @param matrix 2D array storing the result (need to be allocated before being passed).
         * @param n size of the matrix (which is square).
         * @param first_col index of the first column to process.
         * @param n_col number of columns to compute.
         * @param n_mpi_proc number of distributed process, should be 1 when called from a non-MPI code (BEWARE: this is not
         * the number of threads)
         * @param mcpp does the result has to be transposed (this has to be the case for the scalapack linear solve).
         */
        void compute_matrix(Array<double> &matrix, int n, int first_col = 0, int n_col = ALL_COLUMNS, int n_mpi_proc = 1,
                bool transpose = DO_NOT_TRANSPOSE) override;

        void add_var (const char* name, double& var) override
        {   this->call_base_method<const char*,double &>(&System_of_eqs::add_var,name,var); }
        void add_var (const char* name, Tensor& var) override
        {   this->call_base_method<const char*,Tensor &>(&System_of_eqs::add_var,name,var);}
        void add_cst (const char* name, double cst) override
        {   this->call_base_method<const char*,double>(&System_of_eqs::add_cst,name,cst);}
        void add_cst (const char* name, const Tensor& cst) override
        {   this->call_base_method<const char *,Tensor const &>(&System_of_eqs::add_cst,name,cst);}
        void add_def (const char* name) override
        {   this->call_base_method<const char*>(&System_of_eqs::add_def,name);}
        void add_def (int dd, const char* name) override
        {   this->call_base_method<int,const char *>(&System_of_eqs::add_def,dd,name);}
        void add_def_global (const char* name) override
        {   this->call_base_method<const char*>(&System_of_eqs::add_def_global,name);}
        void add_def_global (int dd, const char* name) override
        {   this->call_base_method<int,const char *>(&System_of_eqs::add_def_global,dd,name);}
        void add_ope (const char* name, Term_eq (*pope) (const Term_eq&, Param*), Param* par) override
        {   this->call_base_method<const char *,Term_eq (*)(const Term_eq&,Param *),Param *>(
                &System_of_eqs::add_ope,name,pope,par);}
        void add_ope (const char* name, Term_eq (*pope) (const Term_eq&, const Term_eq&, Param*), Param* par) override
        {   this->call_base_method<const char*,Term_eq (*)(const Term_eq&,const Term_eq&,Param*),Param*>(
                &System_of_eqs::add_ope,name,pope,par);}
        void add_eq_inside (int dom, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,const char *,int,Array<int>**>(
                &System_of_eqs::add_eq_inside,dom,eq,n_cmp,p_cmp);}
        void add_eq_inside (int dom, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,const char *,const List_comp&>(&System_of_eqs::add_eq_inside,dom,eq,list);}
        void add_eq_order (int dom, int order, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,int,const char *,int,Array<int>**>(
                &System_of_eqs::add_eq_order,dom,order,eq,n_cmp,p_cmp);}
        void add_eq_order (int dom, int order, const char* eq,  const List_comp& list) override
        {   this->call_base_method<int,int,const char *,const List_comp&>(
                &System_of_eqs::add_eq_order,dom,order,eq,list);}
        void add_eq_bc (int dom, int bb, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,int,const char*,int,Array<int>**>(
                &System_of_eqs::add_eq_bc,dom,bb,eq,n_cmp,p_cmp);}
        void add_eq_bc (int dom, int bb, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,int,const char*,const List_comp&>(&System_of_eqs::add_eq_bc,dom,bb,eq,list);}
        void add_eq_matching (int dom, int bb, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,int,const char*,int,Array<int>**>(
                &System_of_eqs::add_eq_matching,dom,bb,eq,n_cmp,p_cmp);}
        void add_eq_matching (int dom, int bb, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,int,const char*,const List_comp&>(
                &System_of_eqs::add_eq_matching,dom,bb,eq,list);}
        void add_eq_matching_one_side (int dom, int bb, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,int,const char *,int,Array<int>**>(
                &System_of_eqs::add_eq_matching_one_side,dom,bb,eq,n_cmp,p_cmp);}
        void add_eq_matching_one_side (int dom, int bb, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,int,const char*,const List_comp&>(
                &System_of_eqs::add_eq_matching_one_side,dom,bb,eq,list);}
        void add_eq_matching_non_std (int dom, int bb, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,int,const char*,int,Array<int>**>(
                &System_of_eqs::add_eq_matching_non_std,dom,bb,eq,n_cmp,p_cmp);}
        void add_eq_matching_non_std (int dom, int bb, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,int,const char *,const List_comp&>(
                &System_of_eqs::add_eq_matching_non_std,dom,bb,eq,list);}
        void add_eq_matching_import (int dom, int bb, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,int,const char *,int,Array<int>**>(
                &System_of_eqs::add_eq_matching_import,dom,bb,eq,n_cmp,p_cmp);}
        void add_eq_matching_import (int dom, int bb, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,int,const char *,const List_comp&>(
                &System_of_eqs::add_eq_matching_import,dom,bb,eq,list);}
        void add_eq_full (int dom, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,const char*,int,Array<int>**>(&System_of_eqs::add_eq_full,dom,eq,n_cmp,p_cmp);}
        void add_eq_full (int dom, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,const char*,const List_comp&>(&System_of_eqs::add_eq_full,dom,eq,list);}
        void add_eq_one_side (int dom, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,const char*,int,Array<int>**>(
                &System_of_eqs::add_eq_one_side,dom,eq,n_cmp,p_cmp);}
        void add_eq_one_side (int dom, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,const char*,const List_comp&>(&System_of_eqs::add_eq_one_side,dom,eq,list);}
        void add_eq_matching_exception (int dom, int bb, const char* eq, const Param& par, const char* eq_exception,
                int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,int,const char *,const Param &,const char *,int,Array<int>**>(
                &System_of_eqs::add_eq_matching_exception,dom,bb,eq,par,eq_exception,n_cmp,p_cmp);}
        void add_eq_matching_exception (int dom, int bb, const char* eq, const Param& par, const char* eq_exception,
                const List_comp& list) override
        {   this->call_base_method<int,int,const char*,const Param &,const char *,const List_comp &>(
                &System_of_eqs::add_eq_matching_exception,dom,bb,eq,par,eq_exception,list);}
        void add_eq_order (int dom, const Array<int>& orders, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,const Array<int> &,const char*,int,Array<int>**>(
                &System_of_eqs::add_eq_order,dom,orders,eq,n_cmp,p_cmp);}
        void add_eq_order (int dom, const Array<int>& orders, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,const Array<int>&,const char*,const List_comp&>(
                &System_of_eqs::add_eq_order,dom,orders,eq,list);}
        void add_eq_vel_pot (int dom, int order, const char* eq, const char* const_part) override
        {   this->call_base_method<int,int,const char*,const char*>(
                &System_of_eqs::add_eq_vel_pot,dom,order,eq,const_part);}
        void add_eq_bc (int dom, int bb, const Array<int>& orders, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,int,const Array<int>&,const char*,int,Array<int>**>(
                &System_of_eqs::add_eq_bc,dom,bb,orders,eq,n_cmp,p_cmp);}
        void add_eq_bc (int dom, int bb, const Array<int>& orders, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,int,const Array<int> &,const char*,const List_comp&>(
                &System_of_eqs::add_eq_bc,dom,bb,orders,eq,list);}
        void add_eq_matching (int dom, int bb, const Array<int>& orders, const char* eq, int n_cmp = -1,
                Array<int>** p_cmp=0x0) override
        {   this->call_base_method<int,int,const Array<int>&,const char*,int,Array<int>**>(
                &System_of_eqs::add_eq_matching,dom,bb,orders,eq,n_cmp,p_cmp);}
        void add_eq_matching (int dom, int bb, const Array<int>& orders, const char* eq, const List_comp& list) override
        {   this->call_base_method<int,int,const Array<int> &,const char*,const List_comp&>(
                &System_of_eqs::add_eq_matching,dom,bb,orders,eq,list);}
//        void add_eq_first_integral (int dom, const char* eq, int n_cmp = -1, Array<int>** p_cmp=0x0) override
//        {   this->call_base_method<int,const char*,int,Array<int>**>(
//                &System_of_eqs::add_eq_first_integral,dom,eq,n_cmp,p_cmp);}
        void add_eq_first_integral (int dom_min, int dom_max, const char* integ_part, const char* const_part) override
        {   this->call_base_method<int,int,const char*,const char*>(
                &System_of_eqs::add_eq_first_integral,dom_min,dom_max,integ_part,const_part);}
        void add_eq_mode (int dom, int bb, const char* eq, const Index& pos_cf, double val) override
        {   this->call_base_method<int,int,const char*,const Index &,double>(
                &System_of_eqs::add_eq_mode,dom,bb,eq,pos_cf,val);}
        void add_eq_val_mode (int dom, const char* eq, const Index& pos_cf, double val) override
        {   this->call_base_method<int,const char*,const Index &,double>(
                &System_of_eqs::add_eq_val_mode,dom,eq,pos_cf,val);}
        void add_eq_val (int dom, const char* eq, const Index& pos) override
        {   this->call_base_method<int,const char*,const Index&>(&System_of_eqs::add_eq_val,dom,eq,pos);}
        void add_eq_point (int dom, const char* eq, const Point& MM) override
        {   this->call_base_method<int,const char*,const Point&>(&System_of_eqs::add_eq_point,dom,eq,MM);}

//        Array<double> sec_member() override;
        void compute_nbr_of_conditions() override;
    };

    template<typename... Args>
    System_of_eqs_threaded::System_of_eqs_threaded(std::size_t _max_num_threads,Args && ... args) :
        System_of_eqs{std::forward<Args>(args)...},
        max_num_threads{_max_num_threads},
        thread_specific_systems(max_num_threads-1)
    {
        for(auto & p_syst : thread_specific_systems)
        {
            p_syst.reset(new System_of_eqs{std::forward<Args>(args)...});
        }
    }

    template<typename... Args>
    void System_of_eqs_threaded::call_base_method(void (System_of_eqs::*method)(Args...), Args... args)
    {
        (static_cast<System_of_eqs*>(this)->*static_cast<void (System_of_eqs::*)(Args...)>(method))(args...);
        for(auto & p_syst : thread_specific_systems)
        {
            (static_cast<System_of_eqs&>(*p_syst).*method)(std::forward<Args>(args)...);
        }
    }

}

#endif //__MULTI_THREAD_HPP_
