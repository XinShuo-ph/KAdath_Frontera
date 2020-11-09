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
    //!
    class Solver;

    template<typename T,typename Derived> struct Parameter_base {
        static_assert(!std::is_reference<T>::value,"T cannot be a reference type, consider using pointers");
        using value_type = T;
        using Base = Parameter_base;
        static constexpr value_type get_default_value() {return Derived::default_value;}
        static inline void set(Solver *,Parameter_base const &);
        value_type value;
        constexpr explicit Parameter_base(value_type _value) : value{_value} {}
        constexpr Parameter_base operator=(value_type _value) const {return Parameter_base{_value};}
        constexpr Parameter_base(Parameter_base const &) = default;
        constexpr Parameter_base(Parameter_base&&) noexcept = default;
        constexpr Parameter_base() : value{get_default_value()} {}
    };

    template<typename... Args> struct Setter {static void apply(Solver *,Args && ... ) {};};
    template<typename Arg1,typename... OtherArgs> struct Setter<Arg1,OtherArgs...> {
        static inline void apply(Solver *pthis,Arg1 const & a1,OtherArgs const &... other_args) {
            Arg1::set(pthis,a1);
            Setter<OtherArgs...>::apply(pthis,other_args...);
        }
    };

//    template<> struct Setter<> {static inline void apply(Solver *) {}};


    struct param_output : Parameter_base<bool,param_output> {
        static constexpr bool default_value{true};
        using Base::operator=;
        using Base::Base;
    };
    [[maybe_unused]] constexpr param_output  EnableOutput {true};
    [[maybe_unused]] constexpr param_output  DisableOutput {false};

    struct param_enable_gpu : Parameter_base<bool,param_enable_gpu> {
        static constexpr bool default_value {false};
        using Base::operator=;
        using Base::Base;
    };
    [[maybe_unused]] constexpr param_enable_gpu EnableGPU{true};
    [[maybe_unused]] constexpr param_enable_gpu DisableGPU{false};


    struct param_target_error : Parameter_base<double,param_target_error> {
        static constexpr double default_value{1.e-12};
        using Base::operator=;
        using Base::Base;
    };
    [[maybe_unused]] constexpr param_target_error Tolerance;


    struct param_target_nb_iteration : Parameter_base<int,param_target_nb_iteration> {
        static constexpr int default_value{-1};
        using Base::operator=;
        using Base::Base;
    };
    [[maybe_unused]] constexpr param_target_nb_iteration MaxNbIter;


    struct param_target_duration : Parameter_base<double,param_target_duration> {
        static constexpr double default_value{-1.};
        using Base::operator=;
        using Base::Base;
    };
    [[maybe_unused]] constexpr param_target_duration MaxElapsedTime;


    struct param_target_error_decrease : Parameter_base<double,param_target_error_decrease> {
        static constexpr double default_value {-1.};
        using Base::operator=;
        using Base::Base;
    };
    [[maybe_unused]] constexpr param_target_error_decrease MinImprovement;

    struct param_verbosity : Parameter_base<int,param_verbosity> {
        static constexpr int default_value {1};
        using Base::operator=;
        using Base::Base;
    };
    [[maybe_unused]] constexpr param_verbosity Verbosity;

    struct param_output_stream : Parameter_base<std::ostream*,param_output_stream> {
        static constexpr std::ostream * default_value {&std::cout};
        using Base::operator=;
        using Base::Base;
    };
    [[maybe_unused]] constexpr param_output_stream Output_stream;

    /**
     * Base class for the implementation of functors which solve the system of equations described
     * by an instance of \c System_of_eqs.
     */
    class Solver {
    public:
        using Output_data = System_of_eqs::Output_data;
        enum Stopping_criteria : int {
            fatal_error = -2,
            out_of_memory = -1,
            none = 0,
            tolerance_reached = 1,
            min_error_improvement_reached = 2,
            max_elapsed_time_reached = 3,
            max_nb_iteration_reached = 4
        };

    protected:
        double target_error;
        double current_error;
        double target_error_decrease;
        double current_error_decrease;
        double target_duration;
        double current_duration;
        int target_nb_iteration;
        int current_nb_iteration;
        bool output;
        bool enable_gpu;
        int verbosity;
        std::ostream * output_stream;

        bool master_processus;

    public:
        static constexpr double unknown_value{-1.};
        //! Default constructor.
        Solver() : target_error{param_target_error::default_value}, current_error{std::numeric_limits<double>::max()},
                   target_nb_iteration{param_target_nb_iteration::default_value}, current_nb_iteration{0},
                   target_duration{param_target_duration::default_value}, current_duration{0.},
                   target_error_decrease{param_target_error_decrease::default_value}, current_error_decrease{unknown_value},
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

        Solver & reset_current_values();
        [[nodiscard]] std::pair<bool,Stopping_criteria> check_stopping_criteria() const;

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

        [[nodiscard,maybe_unused]] double get_current_error() const {return current_error;}
        [[nodiscard,maybe_unused]] double get_target_error() const {return target_error;}
        [[nodiscard,maybe_unused]] double get_target_error_decrease() const {return target_error_decrease;}
        [[nodiscard,maybe_unused]] double get_current_error_decrease() const {return current_error_decrease;}
        [[nodiscard,maybe_unused]] double get_target_duration() const {return target_duration;}
        [[nodiscard,maybe_unused]] double get_current_duration() const {return current_duration;}
        [[nodiscard,maybe_unused]] int get_target_nb_iteration() const {return target_nb_iteration;}
        [[nodiscard,maybe_unused]] int get_current_nb_iteration() const {return current_nb_iteration;}
        [[nodiscard,maybe_unused]] bool get_output() const {return output;}
        [[nodiscard,maybe_unused]] bool get_enable_gpu() const {return enable_gpu;}
        [[nodiscard,maybe_unused]] int get_verbosity() const {return verbosity;}
        [[nodiscard,maybe_unused]] std::ostream & get_output_stream() const {return *output_stream;}

        Solver & set_target_error(double new_value) {target_error = new_value; return *this;}
        Solver & set_target_error_decrease(double new_value) {target_error_decrease = new_value; return *this;}
        Solver & set_target_duration(double new_value) {target_duration = new_value; return *this;}
        Solver & set_target_nb_iteration(int new_value) {target_nb_iteration = new_value; return *this;}
        Solver & set_output(bool new_value) {output = new_value; return *this;}
        Solver & set_enable_gpu(bool new_value) {enable_gpu = new_value; return *this;}
        Solver & set_verbosity(int new_value) {verbosity = new_value; return *this;}
        Solver & set_output_stream(std::ostream & new_value) {output_stream = &new_value; return *this;}
    };


    template<> inline void Parameter_base<param_output::value_type,param_output>::set
        (Solver * pthis,Parameter_base const & _param) {pthis->set_output (_param.value);}

    template<> inline void Parameter_base<param_target_error::value_type,param_target_error>::set
        (Solver *pthis,Parameter_base const & _param) {pthis->set_target_error(_param.value);}

    template<> inline void Parameter_base<param_target_duration::value_type,param_target_duration>::set
        (Solver *pthis,Parameter_base const & _param) {pthis->set_target_duration(_param.value);}

    template<> inline void Parameter_base<param_target_error_decrease::value_type,param_target_error_decrease>::set
        (Solver *pthis,Parameter_base const & _param) {pthis->set_target_error_decrease(_param.value);}

    template<> inline void Parameter_base<param_target_nb_iteration::value_type,param_target_nb_iteration>::set
        (Solver *pthis,Parameter_base const & _param) {pthis->set_target_nb_iteration(_param.value);}

    template<> inline void Parameter_base<param_enable_gpu::value_type,param_enable_gpu>::set
        (Solver *pthis,Parameter_base const & _param) {pthis->set_enable_gpu(_param.value);}

    template<> inline void Parameter_base<param_verbosity::value_type,param_verbosity>::set
        (Solver *pthis,Parameter_base const & _param) {pthis->set_verbosity(_param.value);}

    template<> inline void Parameter_base<param_output_stream::value_type,param_output_stream>::set
        (Solver *pthis,Parameter_base const & _param) {pthis->set_output_stream(*_param.value);}

}

#endif //__SOLVERS_HPP_
