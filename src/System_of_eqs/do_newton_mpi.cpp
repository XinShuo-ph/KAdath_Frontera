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
#include "config.h"

#ifdef PAR_VERSION
#include "mpi.h"
#endif

#include "system_of_eqs.hpp"
#include "matrice.hpp"
#include "scalar.hpp"
#include "metric.hpp"
#include "array_math.cpp"

namespace Kadath {

    template<>
    bool System_of_eqs::do_newton<Computational_model::mpi_parallel>(double precision, double& error, std::ostream & os,
            Array<double> */*not used*/ ) {
        bool res;
#ifdef PAR_VERSION
        int bsize  {static_cast<int>(default_block_size)};
        niter++;

        // rank and nproc from MPI :
        int rank, nproc;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        MPI_Comm_size (MPI_COMM_WORLD, &nproc);

        if(rank==0 && niter==1 && display_newton_data)
        {
            display_do_newton_report_header(os,precision);
        }
        //vars_to_terms();
        Array<double> second (sec_member());
        error = max(fabs(second));

        if (error<precision) {
            res = true;
            if(rank==0 && display_newton_data)
            {
                display_do_newton_ending_line(os,precision,error);
                os  << endl;
            }
        }
        else {
            int nn = second.get_size(0);
            if (nn!=nbr_unknowns) {
                cerr << "Number of unknowns is different from the number of equations" << endl;
                cerr << "nbr equations = " << nn << endl;
                cerr << "nbr unknowns  = " << nbr_unknowns << endl;
                abort();
            }

            // Computation in a 1D distributed matrice
            int zero = 0;
            int ictxt_in;
            int nprow_in = 1;
            int npcol_in = nproc;
            sl_init_ (&ictxt_in, &nprow_in, &npcol_in);

            // Get my row and mycol
            int myrow_in, mycol_in;
            blacs_gridinfo_ (&ictxt_in, &nprow_in, &npcol_in, &myrow_in, &mycol_in);
            int nblock_per_proc {(nn/bsize)/nproc};
            int remain_block {(nn/bsize) % nproc};

            while (bsize*nproc>nn)
            {
            }
            if (bsize<1) {
                cerr << "Too many processors in do_newton" << endl;
                abort();
            }
            //adjust block size to avoid idle process when the remaining workload is still significant
            int best_bsize {bsize};
            int best_remain_workload {remain_block};
            while((remain_block <= (nproc/2) && nblock_per_proc<12))
            {
                bsize = div(bsize,2).quot;
                if(bsize != 0) {
                    remain_block = (nn / bsize) % nproc;
                    nblock_per_proc = (nn / bsize) / nproc;
                    if (remain_block > best_remain_workload) {
                        best_remain_workload = remain_block;
                        best_bsize = bsize;
                    }
                }
            }
            //if no suitable block size has been found, the best found is chosen.
            if(bsize<1) bsize = best_bsize;

            int nrowloc_in = numroc_ (&nn, &bsize, &myrow_in, &zero, &nprow_in);
            int ncolloc_in = numroc_ (&nn, &bsize, &mycol_in, &zero, &npcol_in);

            Array<double> matloc_in (ncolloc_in, nrowloc_in);
            int start = bsize*rank;

            Hash_key chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
                                                     nn," | matrix computation ");

            compute_matrix_cyclic(matloc_in,nn,start,bsize,nproc,true);

            // Descriptor of the matrix :
            Array<int> descamat_in(9);
            int info;
            descinit_ (descamat_in.set_data(), &nn, &nn, &bsize, &bsize, &zero, &zero, &ictxt_in, &nrowloc_in, &info);

            // Wait for everybody
            MPI_Barrier(MPI_COMM_WORLD);

            Duration const t_load_matrix {this->stop_chrono(chrono_key)};

            chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
                                            nn," | matrix translation ");
            // Now translate to a 2D cyclic decomposition
            int npcol, nprow;
            split_two_d (nproc, npcol, nprow);
            int ictxt;
            sl_init_ (&ictxt, &nprow, &npcol);

            int myrow, mycol;
            blacs_gridinfo_ (&ictxt, &nprow, &npcol, &myrow, &mycol);

            int nrowloc = numroc_ (&nn, &bsize, &myrow, &zero, &nprow);
            int ncolloc = numroc_ (&nn, &bsize, &mycol, &zero, &npcol);

            Array<double> matloc (ncolloc, nrowloc);
            Array<int> descamat(9);
            descinit_ (descamat.set_data(), &nn, &nn, &bsize, &bsize, &zero, &zero, &ictxt, &nrowloc, &info);

            Cpdgemr2d (nn, nn, matloc_in.set_data(), 1, 1, descamat_in.set_data(), matloc.set_data(), 1, 1, descamat.set_data(), ictxt);

            matloc_in.delete_data();

            // Translate the second member :
            Array<double> secloc (nrowloc);
            for (int row=0 ; row<nn ; row++) {
                int pi = div(row/bsize, nprow).rem;
                int li = int(row/(nprow*bsize));
                int xi = div(row, bsize).rem;
                if ((pi==myrow) && (mycol==0))
                    secloc.set(li*bsize+xi) = second(row);
            }

            // Descriptor of the second member
            Array<int> descsec(9);
            int one = 1;
            descinit_ (descsec.set_data(), &nn, &one, &bsize, &bsize, &zero, &zero, &ictxt, &nrowloc, &info);
            Duration const t_trans_matrix {this->stop_chrono(chrono_key)};

            // Inversion
            Array<int> ipiv (nn);
            chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
                                            nn," | matrix inversion ");
            pdgesv_ (&nn, &one, matloc.set_data(), &one, &one, descamat.set_data(), ipiv.set_data(), secloc.set_data(), &one, &one, descsec.set_data(), &info);
            Duration const t_inv_matrix {this->stop_chrono(chrono_key)};

            chrono_key = this->start_chrono("MPI parallel do newton | problem size = ", nn, " | update ");
            // Get the global solution
            Array<double> auxi(nn);
            auxi = 0.;
            for (int row=0 ; row<nn ; row++) {
                int pi = div(row/bsize, nprow).rem;
                int li = int(row/(nprow*bsize));
                int xi = div(row, bsize).rem;
                if ((pi==myrow) && (mycol==0))
                    auxi.set(row) = secloc(li*bsize+xi);
            }

            Array<double> xx (nn);
            MPI_Allreduce (auxi.set_data(), xx.set_data(), nn, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            blacs_gridexit_ (&ictxt);
            blacs_gridexit_ (&ictxt_in);

            newton_update_vars(xx);

            Duration const t_newton_update {this->stop_chrono(chrono_key)};
            if (rank == 0 && display_newton_data) {
                display_do_newton_iteration(os,
                        {niter,nn,error,t_load_matrix,t_trans_matrix,t_inv_matrix,t_newton_update});
            }
            res = false;
        }
#else
        cerr << "Error : cannot call System_of_eqs::do_newton<mpi_parallel> without MPI. " << endl;
#endif
        return res;
    }


    // Definitions for the linesearch-related methods.

    template<>
    void System_of_eqs::check_positive<Computational_model::mpi_parallel>(double delta)
    {
        if (delta < 0.0)
        {
            cerr << "roundoff error in do_newton_with_linesearch" << endl;
            abort();
        }
    }

    template<>
    void System_of_eqs::check_negative<Computational_model::mpi_parallel>(double delta)
    {
        if (delta >= 0.0)
        {
            cerr << "roundoff error in do_newton_with_linesearch" << endl;
            abort();
        }
    }

    template<>
    double System_of_eqs::compute_f<Computational_model::mpi_parallel>(Array<double> const& second)
    {
        double f(0.0);
        for (int i(0) ; i < second.get_nbr() ; ++i) f += pow(second.get_data()[i],2) ;
        return f;
    }

    template<>
    void System_of_eqs::check_size_VS_unknowns<Computational_model::mpi_parallel>(int n)
    {
        if (n != nbr_unknowns)
        {
            cerr << "Number of unknowns is different from the number of equations" << endl;
            cerr << "nbr equations = " << n << endl;
            cerr << "nbr unknowns  = " << nbr_unknowns << endl;
            abort();
        }
    }

    template<>
    void System_of_eqs::check_bsize<Computational_model::mpi_parallel>(int bsize)
    {
        if (bsize < 1)
        {
            cerr << "Too many processrs in do_newton" << endl;
            abort();
        }
    }

    template<>
    void System_of_eqs::compute_matloc<Computational_model::mpi_parallel>(Array<double>& matloc_in, int nn, int bsize)
    {
#ifdef PAR_VERSION
        int rank, nproc;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        MPI_Comm_size (MPI_COMM_WORLD, &nproc);

        int start(bsize*rank);
        bool done(false);
        int current(0);
        while (!done)
        {
            for (int i(0) ; i < bsize ; ++i)
                if (start+i < nn)
                {
                    Array<double> column(do_col_J(start+i));
                    for (int j(0) ; j < nn ; ++j) matloc_in.set(current,j) = column(j);
                    current++;
                }
            start += nproc*bsize;
            if (start >= nn) done = true;
        }
#else
        cerr << "Error : cannot call System_of_eqs::compute_matloc<mpi_parallel> without MPI." << endl;
#endif
    }

    template<>
    void System_of_eqs::translate_second_member<Computational_model::mpi_parallel>(Array<double>& secloc, Array<double> const& second, int nn, int bsize, int nprow, int myrow, int mycol)
    {
        for (int row(0) ; row < nn ; ++row)
        {
            int pi(div(row/bsize, nprow).rem);
            int li(int(row/(nprow*bsize)));
            int xi(div(row, bsize).rem);
            if ((pi == myrow) and (mycol == 0)) secloc.set(li*bsize+xi) = second(row);
        }
    }

    template<>
    void System_of_eqs::get_global_solution<Computational_model::mpi_parallel>(Array<double>& auxi, Array<double> const& secloc, int nn, int bsize, int nprow, int myrow, int mycol)
    {
        auxi = 0.0;
        for (int row(0) ; row < nn ; ++row)
        {
            int pi(div(row/bsize, nprow).rem);
            int li(int(row/(nprow*bsize)));
            int xi(div(row, bsize).rem);
            if ((pi == myrow) and (mycol == 0)) auxi.set(row) = secloc(li*bsize+xi);
        }
    }

    template<>
    void System_of_eqs::update_fields<Computational_model::mpi_parallel>(double lambda, vector<double> const& old_var_double, vector<Tensor> const& old_var_fields, vector<double> const& p_var_double, vector<Tensor> const& p_var_fields)
    {
#ifdef PAR_VERSION
        int rank;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        if (nvar_double > 0) for (int i(0) ; i < nvar_double ; ++i) *var_double[i] = old_var_double[i] - lambda*p_var_double[i];
        for (int i(0) ; i < nvar ; ++i) *var[i] = old_var_fields[i] - lambda*p_var_fields[i];
#else
        cerr << "Error : cannot call System_of_eqs::update_fields<mpi_parallel> without MPI." << endl;
#endif
    }

    template<>
    void System_of_eqs::compute_old_and_var<Computational_model::mpi_parallel>(Array<double> const& xx, vector<double>& old_var_double, vector<Tensor>& old_var_fields, vector<double>& p_var_double, vector<Tensor>& p_var_fields)
    {
#ifdef PAR_VERSION
        int rank;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        int conte(0);
        espace.xx_to_vars_variable_domains(this, xx, conte);
        if (nvar_double > 0) for (int i(0) ; i < nvar_double ; ++i) old_var_double.push_back(*var_double[i]);
        for (int i(0) ; i < nvar ; ++i) old_var_fields.push_back(Tensor(*var[i]));
        xx_to_vars(xx, conte);
        if (nvar_double > 0) for (int i(0) ; i < nvar_double ; ++i) p_var_double.push_back(*var_double[i]);
        for (int i(0) ; i < nvar ; ++i) p_var_fields.push_back(Tensor(*var[i]));
#else
        cerr << "Error : cannot call System_of_eqs::compute_old_var<mpi_parallel> without MPI." << endl;
#endif
    }

    template<>
    void System_of_eqs::compute_p<Computational_model::mpi_parallel>(Array<double>& xx, Array<double> const& second, int nn)
    {
#ifdef PAR_VERSION
        int bsize(64);
        // rank and nproc from MPI :
        int rank, nproc;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        MPI_Comm_size (MPI_COMM_WORLD, &nproc);

        // Computation in a 1D distributed matrice
        int zero(0);
        int ictxt_in;
        int nprow_in(1);
        int npcol_in(nproc);
        sl_init_(&ictxt_in, &nprow_in, &npcol_in);

        // Get my row and mycol
        int myrow_in, mycol_in;
        blacs_gridinfo_ (&ictxt_in, &nprow_in, &npcol_in, &myrow_in, &mycol_in);
        while (bsize*nproc > nn) bsize = div(bsize,2).quot;
        check_bsize(bsize);
        int nrowloc_in(numroc_(&nn, &bsize, &myrow_in, &zero, &nprow_in));
        int ncolloc_in(numroc_(&nn, &bsize, &mycol_in, &zero, &npcol_in));
        Array<double> matloc_in(ncolloc_in, nrowloc_in);
        compute_matloc(matloc_in, nn, bsize);

        // Descriptor of the matrix :
        Array<int> descamat_in(9);
        int info;
        descinit_(descamat_in.set_data(), &nn, &nn, &bsize, &bsize, &zero, &zero, &ictxt_in, &nrowloc_in, &info);

        // Wait for everybody
        MPI_Barrier(MPI_COMM_WORLD);

        // Now translate to a 2D cyclic decomposition
        int npcol, nprow;
        int ictxt;
        int myrow, mycol;
        split_two_d(nproc, npcol, nprow);
        sl_init_(&ictxt, &nprow, &npcol);
        blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);
        int nrowloc(numroc_(&nn, &bsize, &myrow, &zero, &nprow));
        int ncolloc(numroc_(&nn, &bsize, &mycol, &zero, &npcol));
        Array<double> matloc(ncolloc, nrowloc);
        Array<int> descamat(9);
        descinit_(descamat.set_data(), &nn, &nn, &bsize, &bsize, &zero, &zero, &ictxt, &nrowloc, &info);
        Cpdgemr2d(nn, nn, matloc_in.set_data(), 1, 1, descamat_in.set_data(), matloc.set_data(), 1, 1, descamat.set_data(), ictxt);
        matloc_in.delete_data();

        // Translate the second member :
        Array<double> secloc(nrowloc);
        translate_second_member(secloc, second, nn, bsize, nprow, myrow, mycol);

        // Descriptor of the second member
        Array<int> descsec(9);
        int one(1);
        descinit_(descsec.set_data(), &nn, &one, &bsize, &bsize, &zero, &zero, &ictxt, &nrowloc, &info);

        // Inversion
        Array<int> ipiv(nn);
        pdgesv_(&nn, &one, matloc.set_data(), &one, &one, descamat.set_data(), ipiv.set_data(), secloc.set_data(), &one, &one, descsec.set_data(), &info);

        // Get the global solution
        Array<double> auxi(nn);
        get_global_solution(auxi, secloc, nn, bsize, nprow, myrow, mycol);
        MPI_Allreduce(auxi.set_data(), xx.set_data(), nn, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        blacs_gridexit_(&ictxt);
        blacs_gridexit_(&ictxt_in);
#else
        cerr << "Error : cannot call System_of_eqs::compute_p<mpi_parallel> without MPI." << endl;
#endif
    }

    template<>
    bool System_of_eqs::do_newton_with_linesearch<Computational_model::mpi_parallel>(double precision, double& error,
                                                                         int ntrymax, double stepmax, std::ostream & os)
    {
#ifdef PAR_VERSION
        // Numerical recipes 2007, section 9.7.1
        bool res(false);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        Array<double> second(sec_member());  // rhs of J.x = F
        error = max(fabs(second));
        if (rank == 0) os << "Error init = " << error << endl;
        if (error < precision) res = true;
        else
        {
            int nn(second.get_size(0));
            check_size_VS_unknowns(nn);
            if (rank == 0) os << "Size of the system " << nn << endl;

            Array<double> p(nn);             // solution of J.x = F
            double fold, f, slope, f2(0.0);
            double lambda(1.0), lambdatmp(0.0), lambda2(0.0);
            const double lambdamin(1.0e-10);
            const double alpha(1.0e-4);

            fold = compute_f(second);    // old value of f = 0.5*F.F
            slope = -2.0*fold;
            check_negative(slope);
            compute_p<Computational_model::mpi_parallel>(p, second, nn);    // all the parallel inversion of the jacobian is crammed into this function
            double pmax(max(p));
            if (rank == 0) os << "max(newton step) = " << pmax << endl;

            if (pmax > stepmax) p *= stepmax/pmax;      // rescale in case p is too large
            vector<double> old_var_double, p_var_double;
            vector<Tensor> old_var_fields, p_var_fields;
            compute_old_and_var(p, old_var_double, old_var_fields, p_var_double, p_var_fields);
            for (int itry(1) ; itry <= ntrymax ; ++itry)
            {
                update_fields(lambda, old_var_double, old_var_fields, p_var_double, p_var_fields);    // var now updated
                second = sec_member();
                f = compute_f(second);
                if (f <= 5.0e-13) break; // in double precision, linesearch in spoiled by numerical error as soon as the error is close to 1.0e-7 (f \propto error^2)
                if (f <= fold + alpha*lambda*slope or lambda <= lambdamin) break;
                if (rank == 0) os << "line search " << itry << " ====> ";
                if (lambda == 1.0) lambdatmp = -0.5*slope/(f - fold - slope);
                else
                {
                    double rhs1(f - fold - lambda*slope);
                    double rhs2(f2 - fold - lambda2*slope);
                    double a((rhs1/pow(lambda,2) - rhs2/pow(lambda2,2))/(lambda-lambda2));
                    double b((-lambda2*rhs1/pow(lambda,2) + lambda*rhs2/pow(lambda2,2))/(lambda-lambda2));
                    if (a == 0.0) lambdatmp = -0.5*slope/b;
                    double delta(pow(b,2) - 3.0*a*slope);
                    check_positive<Computational_model::mpi_parallel>(delta);
                    lambdatmp = (-b + sqrt(delta))/(3.0*a);
                    if (lambdatmp > 0.5*lambda) lambdatmp = 0.5*lambda;
                }
                lambda2 = lambda;
                f2 = f;
                lambda = max(lambdatmp, 0.1*lambda);
                if (rank == 0) os << "lambda = " << lambda << endl;
                if (itry == ntrymax) update_fields(lambda, old_var_double, old_var_fields, p_var_double, p_var_fields);
            }
        }
        return res;
#else
        cerr << "Error : cannot call System_of_eqs::do_newton_with_linesearch<mpi_parallel> without MPI." << endl;
        return false;
#endif
    }
}
