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
#include "array_math.hpp"




namespace Kadath {


    template<>
    bool System_of_eqs::do_newton<Computational_model::gpu_mpi_parallel>(double precision, double& error,std::ostream &os,
                                                                       Array<double> * copy_matrix)
    {
        bool res;
#ifdef PAR_VERSION
#ifdef ENABLE_GPU_USE
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
            int const nb_cols_per_proc {nn / nproc};
            int const remaining_cols {nn % nproc};
            int const local_nb_cols {rank < remaining_cols ? nb_cols_per_proc + 1 : nb_cols_per_proc};
            int const local_col_start_idx {rank * nb_cols_per_proc + (rank < remaining_cols ? rank : remaining_cols)};
//	    if(rank==0)
//            {
//                std::cout << "Computing " << nn << 'x' << nn << " matrix with " << nproc << " process." << std::endl;
//                if(remaining_cols==0) std::cout << "Process 0 to " << nproc-1 << " computes " << nb_cols_per_proc << " each." << std::endl;
//	        else std::cout << "- process 0 to " << remaining_cols-1 << " : " << nb_cols_per_proc+1 << " columns " << std::endl
//		    << "- process " << remaining_cols << " to " << nproc-1 << " : " << nb_cols_per_proc << " columns " << std::endl;
//            }
            
            Hash_key chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
                                                     nn," | matrix computation ");
            Array<double> local_matrix_slice (local_nb_cols,nn);
            for(int j{0};j<local_nb_cols;j++)
            {
                int const jj {j+local_col_start_idx};
                Array<double> const colj {do_col_J(jj)};
                for(int i{0};i<nn;i++) local_matrix_slice.set(j, i) = colj(i);
            }
	
    	    std::unique_ptr<Array<double>> full_matrix{nullptr};
	        if(rank==0) full_matrix.reset(new Array<double>(nn,nn));

            std::vector<int> recvcount(nproc,nb_cols_per_proc*nn),displs(nproc,0);
            assert(recvcount.at(rank) == nn*local_nb_cols);
            for(auto i = recvcount.begin();i!=recvcount.begin()+remaining_cols;i++) (*i) += nn;
            for(int i=1;i<nproc;i++) displs.at(i) = displs.at(i-1) + recvcount.at(i-1);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gatherv(local_matrix_slice.get_data(), nn * local_nb_cols, MPI_DOUBLE, (rank==0 ? full_matrix->set_data() : nullptr),
                recvcount.data(),displs.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
            
            Duration const t_load_matrix {this->stop_chrono(chrono_key)};

            Array<double> xx(nn);
            Duration  t_trans_matrix,t_inv_matrix;
            if(rank==0) {
                chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
                                                nn, " | matrix translation ");
                Magma_matrix magma_mat{*full_matrix,nn};
        		full_matrix.reset(nullptr);
                // Translate the second member :
                Magma_array second_member{second};
                t_trans_matrix = this->stop_chrono(chrono_key);

                chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
                                                nn, " | matrix inversion ");
                magma_mat.solve(second_member);
                t_inv_matrix = this->stop_chrono(chrono_key);
                for (int i = 0; i < nn; i++) xx.set(i) = second_member[i];
            }

            chrono_key = this->start_chrono("MPI parallel do newton | problem size = ", nn, " | update ");
            MPI_Bcast(xx.set_data(),nn,MPI_DOUBLE,0,MPI_COMM_WORLD);

            newton_update_vars(xx);

            Duration const t_newton_update{this->stop_chrono(chrono_key)};

            if (rank == 0 && display_newton_data) {
                display_do_newton_iteration(os,
                                            {niter,nn,error,t_load_matrix,t_trans_matrix,t_inv_matrix,t_newton_update});
            }
            res = false;
        }
#endif //#ifdef ENABLE_GPU_USE
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

