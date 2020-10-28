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


    struct output__ {
        static const bool default_value;
        bool value;
        output__ (bool x = default_value) : value{x} {}
        output__ & operator=(bool x) {value = x;return *this;}
    };
    static output__  EnableOutput ;
    static output__  DisableOutput ;


    template<typename... Args> struct Setter {static void apply(Solver *pthis,Args const & ... args);};
    template<typename Arg1,typename... OtherArgs> struct Setter<Arg1,OtherArgs...> {
        static void apply(Solver *pthis,Arg1 const &a1,OtherArgs const &... other_args) {
            Setter<Arg1>::apply(pthis,a1);
            Setter<OtherArgs...>::apply(pthis,other_args...);
        }
    };

    template<> struct Setter< output__ > {
        static void apply(Solver * pthis,output__ const & x) ;
    };

#define decl_algorithm_boolean(member_identifier,keyword)                                                                                 \
    protected:                                                                              \
        bool member_identifier;                                                             \
    public:                                                                                 \
        bool get_ ## member_identifier() const {return member_identifier;}                  \
        Solver & set_ ## member_identifier(bool new_value) {                                \
            member_identifier = new_value;                                                  \
            return *this;                                                                   \
        }                                                                                   \
                                                                                             \

        struct target_ ## member_identifier ##__ {
            using value_type = type;
            static type const default_value;
            type value;
            target_ ## member_identifier ##__ (type x = default_value) : value{x} {}
            target_ ## member_identifier ##__ & operator=(type x) {
                value = x;
                return *this;
            }
        };
        static target_ ## member_identifier ##__ keyword;
        template<> struct Setter< target_ ## member_identifier ##__ > {
            static void apply(Solver *pthis,target_ ## member_identifier ##__ const & x) {
                pthis->set_##target_ ## member_identifier  (x.value);
            }
        };

#define decl_algorithm_parameter(member_identifier,type)                                    \
    protected:                                                                              \
        type target_ ## member_identifier ;                                                 \
        type current_ ## member_identifier ;                                                \
    public:                                                                                 \
        type get_##target_ ## member_identifier  () const {                                 \
            return target_ ## member_identifier ;                                           \
        }                                                                                   \
        Solver & set_##target_ ## member_identifier  (type new_value) {                     \
            target_ ## member_identifier  = new_value;                                      \
            return *this;                                                                   \
        }                                                                                   \
        type get_current_ ## member_identifier () const {                                   \
            return current_ ## member_identifier;                                           \
        }\


#define init_param_setter_buffer(member_identifier,keyword,defvalue)                    \
    Kadath::Solver::target_ ## member_identifier ##__ ::value_type const                \
        Kadath::Solver::target_ ## member_identifier ## __ ::default_value {defvalue};  \
    Kadath::Solver::target_ ## member_identifier ## __ keyword {};         \



    /**
     * Base class for the implementation of functors which solve the system of equations described
     * by an instance of \c System_of_eqs.
     */
    class Solver {
    public:

        decl_algorithm_parameter(error,double);
        decl_algorithm_parameter(nb_iteration,int);
        decl_algorithm_parameter(duration,double);
        decl_algorithm_parameter(error_decrease,double);
        decl_algorithm_boolean(output, Output);

    public:
        static constexpr double unknown_value{-1.};
        //! Default constructor.
        Solver() : target_error{Tolerance.value}, current_error{std::numeric_limits<double>::max()},
                    target_nb_iteration{NbIter.value}, current_nb_iteration{0},
                    target_duration{ElapsedTime.value}, current_duration{0.},
                    target_error_decrease{MinImprovement.value}, current_error_decrease{unknown_value},
                    output{true} {}

        /**
         * Variadic template constructor. Allows to initialize the algorithm parameters using any subset of the
         * available keywords, in any order with the syntax <keyword> = <value>.
         * @tparam Parameter_types catch the parameters type
         * @param parameters may be any subset of the folowing list : Tolerance, NbIter, ElapsedTime, MinImprovement,
         * DisplayOutput, DisableOutput. Except the two last ones, as mentionned above, they may be followed by
         * an assignment operator with the desired value. Exemple :
         * \c Solve(\c NbIter = 10,\c DisabledOutput,\c ElapsedTime = 3600) will build a solver allowing up to 10
         * iterations, without output displayed, and that will stop after one hour of runtime if the max number of
         * iterations or default tolerance are not reached before that.
         */
        template<typename... Parameter_types> Solver(Parameter_types const & ... parameters) : Solver{} {
            this->set(parameters...);
        }

        template<typename... Args> Solver & set(Args const & ... args) {
            Setter<Args...>::apply(this,args...); return *this;
        }

        virtual bool operator()(System_of_eqs & system);

    };

    decl_solver_setting_buffer(error,double,Tolerance);
}

#endif //__SOLVERS_HPP_
