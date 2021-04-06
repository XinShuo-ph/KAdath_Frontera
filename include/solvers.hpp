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
#ifndef __SOLVERS_HPP_
#define __SOLVERS_HPP_

#include <chrono>
#include <numeric>
#include "system_of_eqs.hpp"

namespace Kadath {
    class Solver;

    /**
     * Base class use to declare parameter for the \c Solver class (and its derived one if there is).
     * Beware : this is meant to provide order-free variadic setters, not to be used as data member (see
     * the solver class).
     * @tparam T
     * @tparam Derived
     */
    template<typename T,typename Derived> struct Parameter_base {
        static_assert(!std::is_reference<T>::value,"T cannot be a reference type, consider using pointers");
        using value_type = T;
        using Base = Parameter_base;
        //! Returns the default valoue of the parameter.
        static constexpr value_type get_default_value() {return Derived::default_value;}
        //! Set the parameter value in the appropriate data member of the \c Solver instance pointed to.
        static inline void set(Solver *,Parameter_base const &);
        //! Value of the parameter.
        value_type value;
        //! Constructor.
        constexpr explicit Parameter_base(value_type _value) : value{_value} {}
        //! Assignment operator hack (this is not an assignment !).
        constexpr Parameter_base operator=(value_type _value) const {return Parameter_base{_value};}
        constexpr Parameter_base(Parameter_base const &) = default;
        constexpr Parameter_base(Parameter_base&&) noexcept = default;
        constexpr Parameter_base() : value{get_default_value()} {}
    };

    //! Helper struct used for the partial specialization of a variadic setter.
    template<typename... Args> struct Setter {
        static void apply(Solver *,Args && ... ) {}
    };
    /**
     * Specialization for parameter packs of length greater than or equal to 1.
     * @tparam Arg1 First argument type.
     * @tparam OtherArgs Parameter types of the eventual other arguments.
     */
    template<typename Arg1,typename... OtherArgs> struct Setter<Arg1,OtherArgs...> {
        /**
         * Function that sets the values of argument to the corresponding data members within the \c Solver
         * instance pointed to.
         * @param pthis \Solover instance being parametrized.
         * @param a1 first parameter value
         * @param other_args parameter pack of the the eventual other arguments values.
         */
        static inline void apply(Solver *pthis,Arg1 const & a1,OtherArgs const &... other_args) {
            Arg1::set(pthis,a1);
            Setter<OtherArgs...>::apply(pthis,other_args...);
        }
    };


    //! Parameter controling output during the algorithm execution.
    struct param_output : Parameter_base<bool,param_output> {
        static constexpr bool default_value{true};
        using Base::operator=;
        using Base::Base;
    };
    //! Global variable emulating keywords for the output controling parameter.
    CXX_17_ATTRIBUTES(maybe_unused) constexpr param_output  EnableOutput {true};
    CXX_17_ATTRIBUTES(maybe_unused) constexpr param_output  DisableOutput {false};

    //! Parameter enabling the use of GPU.
    struct param_enable_gpu : Parameter_base<bool,param_enable_gpu> {
        static constexpr bool default_value {false};
        using Base::operator=;
        using Base::Base;
    };
    //! GPU-enabling keyword emulating variable.
    CXX_17_ATTRIBUTES(maybe_unused) constexpr param_enable_gpu EnableGPU{true};
    //! GPU-disabling keyword emulating variable.
    CXX_17_ATTRIBUTES(maybe_unused) constexpr param_enable_gpu DisableGPU{false};

    //! Tolerance/precision.
    struct param_tolerance : Parameter_base<double,param_tolerance> {
        static constexpr double default_value{1.e-12};
        using Base::operator=;
        using Base::Base;
    };
    //! Tolerance keyword emulating variable.
    CXX_17_ATTRIBUTES(maybe_unused) constexpr param_tolerance Tolerance;


    //! Number of iterations.
    struct param_max_nb_iter : Parameter_base<int,param_max_nb_iter> {
        static constexpr int default_value{-1};
        using Base::operator=;
        using Base::Base;
    };
    //! Keyword.
    CXX_17_ATTRIBUTES(maybe_unused) constexpr param_max_nb_iter MaxNbIter;


    struct param_max_elapsed_time : Parameter_base<double,param_max_elapsed_time> {
        static constexpr double default_value{-1.};
        using Base::operator=;
        using Base::Base;
    };
    CXX_17_ATTRIBUTES(maybe_unused) constexpr param_max_elapsed_time MaxElapsedTime;


    struct param_min_improvement : Parameter_base<double,param_min_improvement> {
        static constexpr double default_value {-1.};
        using Base::operator=;
        using Base::Base;
    };
    CXX_17_ATTRIBUTES(maybe_unused) constexpr param_min_improvement MinImprovement;

    struct param_verbosity : Parameter_base<int,param_verbosity> {
        static constexpr int default_value {1};
        using Base::operator=;
        using Base::Base;
    };
    CXX_17_ATTRIBUTES(maybe_unused) constexpr param_verbosity Verbosity;

    struct param_output_stream : Parameter_base<std::ostream*,param_output_stream> {
        static constexpr std::ostream * default_value {&std::cout};
        using Base::operator=;
        using Base::Base;
    };
    CXX_17_ATTRIBUTES(maybe_unused) constexpr param_output_stream Output_stream;

    /**
     * Base class for the implementation of functors which solve the system of equations described
     * by an instance of \c System_of_eqs.
     * \todo it would be better to make a class hierarchy for stopping criteria and simply have a vector of the base instead of all those data members...
     */
    class Solver {
    public:
        using Output_data = System_of_eqs::Output_data;
        //! Result code for the algorithm.
        enum Stopping_criteria : int {
            //! Bad things happened...
            fatal_error = -2,
            //! Algorithm stopped due to memory limitations.
            out_of_memory = -1,
            //! Default value (should not be returned).
            none = 0,
            //! Success with the desired precision reached.
            tolerance_reached = 1,
            //! Algorithm stopped due to too small progression toward error decrease (user setting).
            min_error_improvement_reached = 2,
            //! Algorithm stopped after exceeding the time limit set by the user.
            max_elapsed_time_reached = 3,
            //! Algorithm stopped after exceeding the maximal number of iterations.
            max_nb_iteration_reached = 4
        };

    protected:
        //! Desired precision (i.e. tolerance / error).
        double tolerance;
        //! Currently reached precision/error.
        double current_error;
        //! Minimal error decrease value (disabled if negative).
        double min_improvement;
        //! Current error decrease.
        double current_improvement;
        //! Maximal allowed time.
        double max_elapsed_time;
        //! Current elapsed time.
        double current_elapsed_time;
        //! Maximal allowed number of iterations.
        int max_nb_iter;
        //! Current number of iterations done.
        int current_nb_iter;
        //! Enables/disables the output during the algorithm progression.
        bool output;
        //! Enables/disables GPU use for linear solvers.
        bool enable_gpu;
        //! Verbosity level.
        int verbosity;
        //! STL output stream to which the output is sent.
        std::ostream * output_stream;

        //! Boolean indicating if this proc has null rank or not (mainly to avoid output duplication when using MPI).
        bool master_processus;

    public:
        static constexpr double unknown_value{-1.};
        //! Default constructor.
        Solver() : tolerance{param_tolerance::default_value}, current_error{std::numeric_limits<double>::max()},
                   max_nb_iter{param_max_nb_iter::default_value}, current_nb_iter{0},
                   max_elapsed_time{param_max_elapsed_time::default_value}, current_elapsed_time{0.},
                   min_improvement{param_min_improvement::default_value}, current_improvement{unknown_value},
                   output{param_output::default_value}, enable_gpu{param_enable_gpu::default_value},
                   verbosity{param_verbosity::default_value}, output_stream{param_output_stream::default_value},
                   master_processus{true}
       {
#ifdef PAR_VERSION
           int processus_rank{0};
           MPI_Comm_rank(MPI_COMM_WORLD,&processus_rank);
           master_processus = (processus_rank == 0);
#endif
       }

        /**
         * Variadic template constructor. Allows to initialize the algorithm parameters using any subset of the
         * available keywords, in any order with the syntax <keyword> = <value>.
         * @tparam Parameter_types catch the parameters type
         * @param parameters may be any subset of the folowing list : Tolerance, MaxNbIter, MaxElapsedTime, MinImprovement,
         * DisplayOutput, DisableOutput. Except the two last ones, as mentionned above, they may be followed by
         * an assignment operator with the desired value. Exemple :
         * \c Solve(\c MaxNbIter = 10,\c DisabledOutput,\c MaxElapsedTime = 3600) will build a solver allowing up to 10
         * iterations, without output displayed, and that will stop after one hour of runtime if the max number of
         * iterations or default tolerance are not reached before that.
         */
        template<typename... Parameter_types> explicit Solver(Parameter_types const & ... parameters) : Solver{} {
            this->set(parameters...);
        }

        /**
         * Variadic template setter. Allows to set parameters on an already built object.
         * @tparam Args
         * @param args
         * @return
         */
        template<typename... Args> Solver & set(Args const & ... args) {
            Setter<Args...>::apply(this,args...); return *this;
        }

        /**
         * Check if one of the set stopping criteria is met.
         * @return true if a stopping criteria is triggered, and the code of the stopping criteria.
         */
        Solver & reset_current_values();
        CXX_17_ATTRIBUTES(nodiscard) std::pair<bool,Stopping_criteria> check_stopping_criteria() const;

        virtual Stopping_criteria operator()(System_of_eqs & system);

        /**
         * Sends the header lines for the Newton-Raphson iterations data.
         * @param os stream to send the header to.
         * @param precision is the tolerance stopping criterion for the algorithm.
         */
        void display_do_newton_report_header();

        /**
         * Sends the last lines for the Newton-Raphson iterations data.
         * @param success true if algorithm ended with success, false in the other cases.
         * @param sc stopping criteria flag indicating which criteria stopped the algorithm.
         */
        void display_do_newton_ending_line(Stopping_criteria sc);
        /**
         * Format end sends the Newton-Raphson iteration data.
         * @param os stream to send the data to.
         * @param data data to display.
         */
        void display_do_newton_iteration(const Output_data & data);

        //! Displays informations about the solver (mainly parameters).
        virtual void display(std::ostream & os) const;

        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) double get_current_error() const {return current_error;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) double get_tolerance() const {return tolerance;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) double get_min_improvement() const {return min_improvement;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) double get_current_improvement() const {return current_improvement;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) double get_max_elapsed_time() const {return max_elapsed_time;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) double get_current_elapsed_time() const {return current_elapsed_time;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) int get_max_nb_iter() const {return max_nb_iter;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) int get_current_nb_iter() const {return current_nb_iter;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) bool get_output() const {return output;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) bool get_enable_gpu() const {return enable_gpu;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) int get_verbosity() const {return verbosity;}
        CXX_17_ATTRIBUTES(nodiscard,maybe_unused) std::ostream & get_output_stream() const {return *output_stream;}

        Solver & set_tolerance(double new_value) { tolerance = new_value; return *this;}
        Solver & set_min_improvement(double new_value) { min_improvement = new_value; return *this;}
        Solver & set_max_elapsed_time(double new_value) { max_elapsed_time = new_value; return *this;}
        Solver & set_max_nb_iter(int new_value) { max_nb_iter = new_value; return *this;}
        Solver & set_output(bool new_value) {output = new_value; return *this;}
        Solver & set_enable_gpu(bool new_value) {enable_gpu = new_value; return *this;}
        Solver & set_verbosity(int new_value) {verbosity = new_value; return *this;}
        Solver & set_output_stream(std::ostream & new_value) {output_stream = &new_value; return *this;}
    };


    template<> inline void Parameter_base<param_output::value_type,param_output>::set
        (Solver * pthis,Parameter_base const & _param) {pthis->set_output (_param.value);}

    template<> inline void Parameter_base<param_tolerance::value_type,param_tolerance>::set
        (Solver *pthis,Parameter_base const & _param) { pthis->set_tolerance(_param.value);}

    template<> inline void Parameter_base<param_max_elapsed_time::value_type,param_max_elapsed_time>::set
        (Solver *pthis,Parameter_base const & _param) { pthis->set_max_elapsed_time(_param.value);}

    template<> inline void Parameter_base<param_min_improvement::value_type,param_min_improvement>::set
        (Solver *pthis,Parameter_base const & _param) { pthis->set_min_improvement(_param.value);}

    template<> inline void Parameter_base<param_max_nb_iter::value_type,param_max_nb_iter>::set
        (Solver *pthis,Parameter_base const & _param) { pthis->set_max_nb_iter(_param.value);}

    template<> inline void Parameter_base<param_enable_gpu::value_type,param_enable_gpu>::set
        (Solver *pthis,Parameter_base const & _param) {pthis->set_enable_gpu(_param.value);}

    template<> inline void Parameter_base<param_verbosity::value_type,param_verbosity>::set
        (Solver *pthis,Parameter_base const & _param) {pthis->set_verbosity(_param.value);}

    template<> inline void Parameter_base<param_output_stream::value_type,param_output_stream>::set
        (Solver *pthis,Parameter_base const & _param) {pthis->set_output_stream(*_param.value);}

}

#endif //__SOLVERS_HPP_
