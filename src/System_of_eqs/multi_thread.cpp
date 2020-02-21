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

#include "multi_thread.hpp"


//Kadath::Array<double> Kadath::System_of_eqs_threaded::sec_member()
//{
//    for(auto & p_syst : thread_specific_systems)
//    {
//        p_syst->sec_member();
//    }
//    return System_of_eqs::sec_member();
//}

void Kadath::System_of_eqs_threaded::compute_matrix(Array<double> &matrix, int n, int first_col, int n_col,
        int n_mpi_proc, bool transpose)
{
    if(n_mpi_proc>1)
    {
        std::cerr << "Error : System_of_eqs_threaded::compute_matrix(&matrix="<< &matrix << ", n=" << n
                  << ", first_col=" << first_col << ", n_col="<< n_col << ", n_mpi_proc=" << n_mpi_proc <<", transpose="
                  << (transpose ? "true" : "false") << " at line " << __LINE__ << " in file " << __FILE__ << std::endl;
        std::cerr << " ==> mixed MPI/multithread parallelism is not handled yet." << std::endl;
        throw runtime_error{"Attempt to mix MPI with multi-threading."};
    }
    assert(first_col==0);
    assert(n_col ==ALL_COLUMNS);
    std::vector<System_of_eqs*> ptr_to_thread_soe = {this};
    for(auto & soe : thread_specific_systems) ptr_to_thread_soe.push_back(&(*soe));
    assert(ptr_to_thread_soe.size() == max_num_threads);
    n_col = (n_col == ALL_COLUMNS ? n : n_col);
    std::size_t block_size{default_block_size};
    while(block_size * max_num_threads > n_col) block_size >>= 1;
    if( block_size == 0)
    {
        std::cerr << "Error : System_of_eqs_threaded::compute_matrix(&matrix="<< &matrix << ", n=" << n
            << ", first_col=" << first_col << ", n_col="<< n_col << ", n_mpi_proc=" << n_mpi_proc <<", transpose="
            << (transpose ? "true" : "false") << " at line " << __LINE__ << " in file " << __FILE__ << std::endl;
        std::cerr << " ==> to many threads (max_num_threads=" << max_num_threads << " compared to the number of columns"
                    " to process (n_col= " << n_col << ")." << std::endl;
        throw std::runtime_error{"Too many threads."};
    }
    std::deque<std::thread> threads{};
    for(std::size_t k {0}; k<max_num_threads; k++)
    {
        System_of_eqs * p2soe {ptr_to_thread_soe.at(k)};
        threads.emplace_back([this,p2soe,&matrix,n,k,block_size,transpose](){p2soe->System_of_eqs::compute_matrix(
                matrix, n,k*block_size,block_size,max_num_threads,transpose);});
//        threads.emplace_back(&System_of_eqs::compute_matrix,p2soe,matrix,n,static_cast<int>(k*block_size),
//                static_cast<int>(block_size),static_cast<int>(max_num_threads),transpose);
    }
    for(auto & thread : threads) thread.join();
}

void Kadath::System_of_eqs_threaded::compute_nbr_of_conditions()
{
    System_of_eqs::compute_nbr_of_conditions();
    for(auto & p_syst : thread_specific_systems)
    {
        p_syst->compute_nbr_of_conditions();
    }
}
