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

#include "solvers.hpp"


namespace Kadath {

    Solver::Stopping_criteria Solver::operator()(System_of_eqs &system) {
        reset_current_values();
        if(target_nb_iteration < 0 && target_error < 0. && target_duration < 0. && target_error_decrease < 0.) {
            if(system.get_mpi_proc_rank() == 0 && verbosity>0)
                std::cerr << "WARNING: No stopping criteria set, high risk of endless loop !" << std::endl;
        }
        std::pair<bool,Stopping_criteria> algo_state {false,none};
        auto start_time = std::chrono::steady_clock::now();
        bool converged {false};
        try{
            while(!algo_state.first && !converged) {
                double previous_error {current_error};
#ifdef PAR_VERSION
#ifdef ENABLE_GPU_USE
                if(enable_gpu) {
                    converged = system.do_newton<Computational_model::gpu_mpi_parallel>(target_error >= 0. ? target_error : 0.,
                                                                                        current_error,output);
                }
                else {
                    converged = system.do_newton<Computational_model::mpi_parallel>(target_error >= 0. ? target_error : 0.,
                                                                                    current_error, output);
                }
#else //ifdef ENABLE_GPU_USE
                converged = system.do_newton<Computational_model::mpi_parallel>(target_error >= 0. ? target_error : 0.,
                                                                                current_error,output);
#endif //ifdef ENABLE_GPU_USE
#else  //ifdef PAR_VERSION
                converged = system.do_newton<Computational_model::sequential>(target_error >= 0. ? target_error : 0.,
                                                                                current_error,output && verbosity>0);
#endif //ifdef PAR_VERSION
                current_nb_iteration++;
                current_error_decrease = previous_error - current_error ;
                std::chrono::steady_clock::duration elapsed {std::chrono::steady_clock::now() - start_time};
                current_duration = std::chrono::duration_cast<std::chrono::duration<double>>(elapsed).count();
                algo_state = check_stopping_criteria();
            }
        }
        catch (std::bad_alloc & e) {
            if(verbosity)
                std::cerr << "ERROR: Iterating over do_newton() lead to the following exception :\n"
                          << e.what() << std::endl;
            return out_of_memory;
        }
        catch(std::exception & e) {
            if(verbosity)
                std::cerr << "ERROR: Iterating over do_newton() lead to the following exception :\n"
                          << e.what() << std::endl;
            return fatal_error;
        }
        system.finalize_profiling();
        if(verbosity > 1) {
            profiling_report(system,system.get_output_stream());
        }
#ifdef ALL_CHECKS_ENABLED
        if(algo_state.second == tolerance_reached && !converged) {
            if(system.get_mpi_proc_rank()==0) {
                if(verbosity)
                    std::cerr << "ERROR: System_of_eqs::do_newton and Solver returned non-matching "
                                 "convergence results.";
                return fatal_error;
            }
        }
#endif
        return algo_state.second;
    }

    std::pair<bool,Solver::Stopping_criteria> Solver::check_stopping_criteria() const {
        if(target_error > 0. && current_error <= target_error) {
            return std::make_pair(true, tolerance_reached);
        }
        else if(target_nb_iteration>=0 && current_nb_iteration >= target_nb_iteration) {
            return std::make_pair(true, max_nb_iteration_reached);
        }
        else if(target_duration >= 0. && current_duration >= target_duration) {
            return std::make_pair(true, max_elapsed_time_reached);
        }
        else if(target_error_decrease >= 0. && current_error_decrease <= target_error_decrease) {
            return std::make_pair(true, min_error_improvement_reached);
        }
        else return std::make_pair(false, none);
    }


    Solver & Solver::reset_current_values() {
        current_error = std::numeric_limits<double>::max();
        current_nb_iteration = 0;
        current_duration = 0.;
        current_error_decrease = -1.;
        return *this;
    }


    void Solver::display(std::ostream &os) const {
        auto EorD = [](bool condition) -> auto {return condition ? "enabled  " : "disabled ";};
        os  << "Kadath::Solver settings :\n";
        os  << "\t - Data output display.................... : " << EorD(output) << std::endl;
        os  << "\t - Verbosity level........................ : " << verbosity << std::endl;
        os  << "\t - Error value-based stop. criteria....... : ";
        if(target_error <= 0.) os << "disabled \n";
        else os << "enabled  -  " << current_error << " / " << target_error << "\n";
        os  << "\t - Error improvement-based stop. criteria. : ";
        if(target_error_decrease <= 0.) os << "disabled \n";
        else os << "enabled  -  " << current_error_decrease << " / " << target_error_decrease << "\n";
        os  << "\t - Iteration-based stop. criteria......... : ";
        if(target_nb_iteration < 0) os << "disabled \n";
        else os << "enabled  -  " << current_nb_iteration << " / " << target_nb_iteration << "\n";
        os  << "\t - Time-based stop. criteria.............. : ";
        if(target_duration <= 0.) os << "disabled \n";
        else os << "enabled  -  " << current_duration << " s / " << target_duration << " s\n";
        os  << "\t - Use GPU for linear solve (if available) : " << EorD(enable_gpu) << std::endl;
#ifndef ENABLE_GPU_USE
        if(enable_gpu) os << "\t WARNING : this version of Kadath has not been build with GPU capabilities,\n"
                          << "\t           this option will be ignored." << std::endl;
#endif
        os << std::endl;
    }
}