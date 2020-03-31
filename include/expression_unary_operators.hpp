//
// Created by sauliac on 27/03/2020.
//

#include <expression.hpp>

#ifndef __EXPRESSION_OPERATORS_HPP_
#define __EXPRESSION_OPERATORS_HPP_

namespace Kadath {
    namespace internal {
        /**
         * Expression type for unary operators or functions of a unique variable.
         * @tparam Operator helper class containing the actual implementation of the operator or function and
         * type-related informations.
         * @tparam Argument_expression the sub-expression the operator or function is applied to.
         */
        template<typename Operator,typename Argument_expression> struct Unary_operator
            : public Expression<Unary_operator<Operator,Argument_expression>,typename Operator::Return_type>
        {
            //! Return value type of this expression.
            using Value_type = typename Operator::Return_type;
            //! Argument value type of the operator.
            using Arg_value_type = typename Operator::Argument_type;
            //! Actual expression type of the argument expression.
            using Arg_expression_type = typename Argument_expression::Expression_type;
            static_assert(std::is_same<Arg_value_type,typename Argument_expression::Value_type>::value);

            //! The sub-expression tree.
            Arg_expression_type const argument;

            //! Constructor.
            explicit Unary_operator(Argument_expression const & _argument) : argument{_argument.cast()} {}
            //! Evaluation of the \c i-th component.
            Value_type evaluate_(int i) const {return Operator::evaluate(argument.evaluate(i));}

            //! Implementation of the \c check_dimensions method, checks the dimensions in the argument.
            bool check_dimensions_() const {return argument.check_dimensions();}
            //! Implementation of the \c get_size method.
            int get_size_() const {return argument.get_size();}
            //! Implementation of the \c get_dimensions method for the array expression case.
            Dim_array get_dimensions_() const {return argument.get_dimensions();}
        };

        /*
         * Unary operators and one variable functions helper structs :
         */

        /**
         * Base class for unary operators and single variable functions helper struct for the instantiation of the
         * \c Unary_operator expression template struct.
         * @tparam A argument type.
         * @tparam R return type.
         */
        template<typename A,typename R=A> struct Unary_operator_helper_base {
            //! Argument type alias.
            using Argument_type = A;
            //! Return type alias.
            using Return_type = R;
        };

        /**
         * Expression template utility helper class for the unary + operator.
         * @tparam A argument type
         */
        template<typename T> struct Unary_add : Unary_operator_helper_base<T> {
            static constexpr T evaluate(T const x) {return + x;}
        };

        /**
         * Expression template utility helper class for the unary - operator.
         * @tparam T argument type
         */
        template<typename T> struct Unary_sub : Unary_operator_helper_base<T> {
            static constexpr T evaluate(T const x) {return - x;}
        };

        /**
         * Expression template utility herpler class for the exponential function (given as example, the other STL
         * functions will be generated with the help of preprocessor macros).
         * @tparam T argument type
         */
        template<typename T> struct Func_exp : Unary_operator_helper_base<T> {
            static inline T evaluate(T const x) {return std::exp(x);}
        };


/**
 * This macro declares and defines, for a given function, the overloads necessary for its use within
 * the expression template optimization.
 */
#define EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(function) \
        template<typename T> struct Func_##function : Unary_operator_helper_base<T> { \
            static inline T evaluate(T const x) {return std:: function (x);} \
        }; \
        template<typename T> Unary_operator<Func_##function <T>,Vector_expression<T>> function (Array<T> const & vector) { \
            return Unary_operator<Func_##function <T>,Vector_expression<T>>{Vector_expression<T>{vector}}; \
        } \
        template<typename Expr> Unary_operator<Func_##function <typename Expr::Value_type>,Expression<Expr>> \
        function(Expression<Expr> const & argument) { \
            return Unary_operator<Func_##function <typename Expr::Value_type>,Expression<Expr>>{argument}; \
        }

        namespace operators{
            //! Unary operator + overload for the \c Array class.
            template<typename T> Unary_operator<Unary_add<T>,Vector_expression<T>> operator+(Array<T> const & vector) {
                return Unary_operator<Unary_add<T>,Vector_expression<T>>{Vector_expression<T>{vector}};
            }
            //! Unary operator + overload for all possible expression types.
            template<typename Expr> Unary_operator<Unary_add<typename Expr::Value_type>,Expression<Expr>>
            operator+(Expression<Expr> const & argument) {
                return Unary_operator<Unary_add<typename Expr::Value_type>,Expression<Expr>>{argument};
            }

            //! Unary operator - overload for the \c Array class.
            template<typename T> Unary_operator<Unary_sub<T>,Vector_expression<T>> operator-(Array<T> const & vector) {
                return Unary_operator<Unary_sub<T>,Vector_expression<T>>{Vector_expression<T>{vector}};
            }
            //! Unary operator - overload for all possible expression types.
            template<typename Expr> Unary_operator<Unary_sub<typename Expr::Value_type>,Expression<Expr>>
            operator-(Expression<Expr> const & argument) {
                return Unary_operator<Unary_sub<typename Expr::Value_type>,Expression<Expr>>{argument};
            }

            //! exp function overload for the \c Array class.
            template<typename T> Unary_operator<Func_exp<T>,Vector_expression<T>> exp(Array<T> const & vector) {
                return Unary_operator<Func_exp<T>,Vector_expression<T>>{Vector_expression<T>{vector}};
            }
            //! exp function overload for all possible expression types.
            template<typename Expr> Unary_operator<Func_exp<typename Expr::Value_type>,Expression<Expr>>
            exp(Expression<Expr> const & argument) {
                return Unary_operator<Func_exp<typename Expr::Value_type>,Expression<Expr>>{argument};
            }
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(sqrt);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(cbrt);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(cos);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(sin);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(tan);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(acos);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(asin);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(atan);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(exp2);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(log);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(log2);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(log1p);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(log10);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(cosh);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(sinh);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(tanh);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(acosh);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(asinh);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(atanh);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(abs);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(fabs);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(erf);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(erfc);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(tgamma);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(lgamma);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(floor);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(ceil);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(trunc);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(round);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(nearbyint);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(rint);
            EXPRESSION_TEMPLATE_FUNCTION_OVERLOADS(logb);
        } // namespace operators

    }// namespace internal

}// namespace Kadath

#endif //__EXPRESSION_OPERATORS_HPP_
