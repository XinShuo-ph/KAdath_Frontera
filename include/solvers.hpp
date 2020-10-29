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
        using value_type = T;
        using Base = Parameter_base;
        static constexpr value_type get_default_value() {return Derived::default_value;}
        static inline void set(Solver *,Parameter_base const &);
        value_type value;
        constexpr Parameter_base(value_type _value = get_default_value()) : value{_value} {}
        constexpr Parameter_base operator=(value_type _value) const {return _value;}
        constexpr Parameter_base(Parameter_base const &) = delete;
        constexpr Parameter_base(Parameter_base&&) = delete;
    };

    template<typename... Args> struct Setter {static void apply(Solver *pthis,Args && ... args);};
    template<typename Arg1,typename... OtherArgs> struct Setter<Arg1,OtherArgs...> {
        static inline void apply(Solver *pthis,Arg1 const & a1,OtherArgs const &... other_args) {
            Arg1::set(pthis,a1);
            Setter<OtherArgs...>::apply(pthis,other_args...);
        }
    };

    template<> struct Setter<> {static inline void apply(Solver *) {}};


    struct param_output : Parameter_base<bool,param_output> {
        static constexpr bool default_value{true};
        using Base::operator=;
    };
    constexpr param_output  EnableOutput {true};
    constexpr param_output  DisableOutput {false};

    struct param_enable_gpu : Parameter_base<bool,param_enable_gpu> {
        static constexpr bool default_value {false};
        using Base::operator=;
    };
    constexpr param_enable_gpu EnableGPU{true};
    constexpr param_enable_gpu DisableGPU{false};


    struct param_target_error : Parameter_base<double,param_target_error> {
        static constexpr double default_value{1.e-12};
        using Base::operator=;
    };
    constexpr param_target_error Tolerance;


    struct param_target_nb_iteration : Parameter_base<int,param_target_nb_iteration> {
        static constexpr int default_value{-1};
        using Base::operator=;
    };
    constexpr param_target_nb_iteration MaxNbIter;


    struct param_target_duration : Parameter_base<double,param_target_duration> {
        static constexpr double default_value{-1.};
        using Base::operator=;
    };
    constexpr param_target_duration MaxElapsedTime;


    struct param_target_error_decrease : Parameter_base<double,param_target_error_decrease> {
        static constexpr double default_value {-1.};
        using Base::operator=;
    };
    constexpr param_target_error_decrease MinImprovement;



    /**
     * Base class for the implementation of functors which solve the system of equations described
     * by an instance of \c System_of_eqs.
     */
    class Solver {
    public:
        enum Stopping_criteria : int {
            none = -1,
            tolerance = 0,
            min_error_improvement = 1,
            max_elapsed_time = 2,
            max_nb_iteration = 3
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

    public:
        static constexpr double unknown_value{-1.};
        //! Default constructor.
        Solver() : target_error{Tolerance.default_value}, current_error{std::numeric_limits<double>::max()},
                   target_nb_iteration{MaxNbIter.default_value}, current_nb_iteration{0},
                   target_duration{MaxElapsedTime.default_value}, current_duration{0.},
                   target_error_decrease{MinImprovement.default_value}, current_error_decrease{unknown_value},
                   output{param_output::default_value}, enable_gpu{param_enable_gpu::default_value} {}

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
        template<typename... Parameter_types> Solver(Parameter_types const & ... parameters) : Solver{} {
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
        std::pair<bool,Stopping_criteria> check_stopping_criteria() const;

        virtual bool operator()(System_of_eqs & system);

        //! Displays informations about the solver (mainly parameters).
        virtual void display(std::ostream & os) const;

        double get_target_error() const {return target_error;}
        double get_current_error() const {return current_error;}
        double get_target_error_decrease() const {return target_error_decrease;}
        double get_current_error_decrease() const {return current_error_decrease;}
        double get_target_duration() const {return target_duration;}
        double get_current_duration() const {return current_duration;}
        int get_target_nb_iteration() const {return target_nb_iteration;}
        int get_current_nb_iteration() const {return current_nb_iteration;}
        bool get_output() const {return output;}
        bool get_enable_gpu() const {return enable_gpu;}

        Solver & set_target_error(double new_value) {target_error = new_value; return *this;}
        Solver & set_target_error_decrease(double new_value) {target_error_decrease = new_value; return *this;}
        Solver & set_target_duration(double new_value) {target_duration = new_value; return *this;}
        Solver & set_target_nb_iteration(int new_value) {target_nb_iteration = new_value; return *this;}
        Solver & set_output(bool new_value) {output = new_value; return *this;}
        Solver & set_enable_gpu(bool new_value) {enable_gpu = new_value; return *this;}

    };


    template<> inline void Parameter_base<bool,param_output>::set(Solver * pthis, Parameter_base const & _param) {
        pthis->set_output (_param.value);
    }

    template<> inline void Parameter_base<double,param_target_error>::set(Solver *pthis, Parameter_base const & _param) {
        pthis->set_target_error(_param.value);
    }

    template<> inline void Parameter_base<double,param_target_duration>::set(Solver *pthis, Parameter_base const & _param) {
        pthis->set_target_duration(_param.value);
    }

    template<> inline void Parameter_base<double,param_target_error_decrease>::set(Solver *pthis, Parameter_base const & _param) {
        pthis->set_target_error_decrease(_param.value);
    }

    template<> inline void Parameter_base<int,param_target_nb_iteration>::set(Solver *pthis, Parameter_base const & _param) {
        pthis->set_target_nb_iteration(_param.value);
    }

    template<> inline void Parameter_base<bool,param_enable_gpu>::set(Solver *pthis,Parameter_base const & _param) {
        pthis->set_enable_gpu(_param.value);
    }

}

#endif //__SOLVERS_HPP_
