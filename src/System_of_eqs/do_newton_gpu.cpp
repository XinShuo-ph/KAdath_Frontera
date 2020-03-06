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

#include "magma_interface.hpp"
#include "system_of_eqs.hpp"
#include "matrice.hpp"
#include "array_math.cpp"




namespace Kadath {


    template<>
    bool System_of_eqs::do_newton<Computational_model::gpu_mpi_parallel>(double precision, double& error,std::ostream &os,
                                                                       Array<double> * copy_matrix)
    {
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

            while (bsize*nproc>nn) {
                bsize = div(bsize,2).quot;
            }
            if (bsize<1) {
                cerr << "Too many processors in do_newton" << endl;
                abort();
            }

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

            // Now translate the matrix back to the master process
            int npcol{1}, nprow{1};
            int ictxt;
            sl_init_ (&ictxt, &nprow, &npcol);

            int myrow, mycol;
            blacs_gridinfo_ (&ictxt, &nprow, &npcol, &myrow, &mycol);

            int nrowloc = numroc_ (&nn, &bsize, &myrow, &zero, &nprow);
            int ncolloc = numroc_ (&nn, &bsize, &mycol, &zero, &npcol);
            if(rank == 0)
            {
                assert(nrowloc == nn && ncolloc == nn);
            } else{
                assert(nrowloc==0 && ncolloc==0);
            }

            Array<double> xx(nn);
            Duration  t_trans_matrix,t_inv_matrix;
            if(rank==0) {
                chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
                                                nn, " | matrix translation ");
                Magma_matrix magma_mat{nn};
                {
                    Array<double> matloc(ncolloc, nrowloc);
                    Array<int> descamat(9);
                    descinit_(descamat.set_data(), &nn, &nn, &bsize, &bsize, &zero, &zero, &ictxt, &nrowloc, &info);

                    Cpdgemr2d(nn, nn, matloc_in.set_data(), 1, 1, descamat_in.set_data(), matloc.set_data(), 1, 1,
                              descamat.set_data(), ictxt);
                    magma_mat = matloc;
                    matloc_in.delete_data();
                }
                // Translate the second member :
                Magma_array second_member{second};
                t_trans_matrix = this->stop_chrono(chrono_key);

                chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
                                                nn, " | matrix inversion ");
                magma_mat.solve(second_member);
                t_inv_matrix = this->stop_chrono(chrono_key);
                for (int i = 0; i < nn; i++) xx.set_data()[i] = second_member[i];
            }

            chrono_key = this->start_chrono("MPI parallel do newton | problem size = ", nn, " | update ");
            MPI_Bcast(xx.set_data(),nn,MPI_DOUBLE,0,MPI_COMM_WORLD);

            blacs_gridexit_(&ictxt);
            blacs_gridexit_(&ictxt_in);

            newton_update_vars(xx);

            Duration const t_newton_update{this->stop_chrono(chrono_key)};

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

    template<>
    bool System_of_eqs::do_newton<Computational_model::gpu_sequential>(double precision, double& error,std::ostream &os,
                                                                   Array<double> * copy_matrix)
    {
#ifdef PAR_VERSION
        int rank;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        if(rank==0) {
#endif
            std::cerr << "Not implemented yet." << std::endl;
            return false;
#ifdef PAR_VERSION
        }
        else {
            return false;
        }
#endif
    }
//#endif


}

