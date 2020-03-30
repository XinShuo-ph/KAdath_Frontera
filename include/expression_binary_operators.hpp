//
// Created by sauliac on 27/03/2020.
//

#include "expression.hpp"

#ifndef __EXPRESSION_BINARY_OPERATORS_HPP_
#define __EXPRESSION_BINARY_OPERATORS_HPP_

namespace Kadath {
    namespace internal {

        /**
         * Expression type for binary operators or functions of two variables.
         * @tparam Operator helper utility class containing the actual call to the operator and type-related infos.
         * @tparam Left_expression left/first argument/variable type.
         * @tparam Right_expression right/second argument/variable type
         */
        template<typename Operator, typename Left_expression, typename Right_expression>
        struct Binary_operator :
                public Expression<Binary_operator<Operator, Left_expression, Right_expression>,
                        typename Operator::Return_type> {
            //! Return value type of this expression.
            using Value_type = typename Operator::Return_type;
            //! Left argument value type.
            using Left_arg_value_type = typename Operator::Left_arg_type;
            //! Left expression actual type.
            using Left_expression_type = typename Left_expression::Expression_type;
            //! Right argument value type.
            using Right_arg_value_type = typename Operator::Right_arg_type;
            //! Right expression actual type.
            using Right_expression_type = typename Right_expression::Expression_type;

            static_assert(std::is_same<Left_arg_value_type, typename Left_expression::Value_type>::value);
            static_assert(std::is_same<Right_arg_value_type, typename Right_expression::Value_type>::value);

            //! Left expression sub-expression tree.
            Left_expression_type const left;
            //! Right expression sub-expression tree.
            Right_expression_type const right;

            //! Constructor.
            explicit Binary_operator(Left_expression const &_left, Right_expression const &_right) :
                    left{_left.cast()}, right{_right.cast()} {}

            //! Implementation of the \c check_dimensions method, checks the dimensions in the arguments.
            bool check_dimensions_() const {
                return left.check_dimensions() && right.check_dimensions() &&
                       left.get_dimensions() == right.get_dimensions();
            }

            //! Implementation of the \c get_size method.
            int get_size_() const { return left.get_size(); }

            //! Implementation of the \c get_dimensions method for the array expression case.
            Dim_array get_dimensions_() const { return left.get_dimensions(); }
        };

        /**
         * Expression type for binary operators or functions of two variables with a scalar value on the left
         * hand side argument / first variable.
         * @tparam Operator helper utility class containing the actual call to the operator and type-related infos.
         * @tparam Right_expression right/second argument/variable type
         */
        template<typename Operator, typename Right_expression>
        struct Binary_operator_left_scalar :
                Expression<Binary_operator_left_scalar<Operator, Right_expression>, typename Operator::Return_type> {
            //! Return value type of this expression.
            using Value_type = typename Operator::Return_type;
            //! Left argument value type.
            using Left_arg_value_type = typename Operator::Left_arg_type;
            //! Right argument value type.
            using Right_arg_value_type = typename Operator::Right_arg_type;
            //! Right expression actual type.
            using Right_expression_type = typename Right_expression::Expression_type;

            static_assert(std::is_same<Right_arg_value_type, typename Right_expression::Value_type>::value);

            //! Left scalar value.
            Left_arg_value_type const left_scalar;
            //! Right expression sub-expression tree.
            Right_expression_type const right_expression;

            //! Constructor.
            explicit Binary_operator_left_scalar(Left_arg_value_type const _left_scalar,
                                                 Right_expression const &_right_expression) :
                    left_scalar{_left_scalar}, right_expression{_right_expression.cast()} {}

            //! Implementation of the \c check_dimensions method, checks the dimensions in the arguments.
            bool check_dimensions_() const { return right_expression.check_dimensions(); }

            //! Implementation of the \c get_size method.
            int get_size_() const { return right_expression.get_size(); }

            //! Implementation of the \c get_dimensions method for the array expression case.
            Dim_array get_dimensions_() const { return right_expression.get_dimensions(); }
        };

        /**
         * Expression type for binary operators or functions of two variables with a scalar value on the right
         * hand side argument / second variable.
         * @tparam Operator helper utility class containing the actual call to the operator and type-related infos.
         * @tparam Left_expression left/first argument/variable type
         */
        template<typename Operator, typename Left_expression>
        struct Binary_operator_right_scalar :
                Expression<Binary_operator_right_scalar<Operator, Left_expression>, typename Operator::Return_type> {
            //! Return value type of this expression.
            using Value_type = typename Operator::Return_type;
            //! Left argument value type.
            using Left_arg_value_type = typename Operator::Left_arg_type;
            //! Right argument value type.
            using Right_arg_value_type = typename Operator::Right_arg_type;
            //! Right expression actual type.
            using Left_expression_type = typename Left_expression::Expression_type;

            static_assert(std::is_same<Right_arg_value_type, typename Left_expression::Value_type>::value);

            //! Right expression sub-expression tree.
            Left_expression_type const left_expression;
            //! Right scalar value.
            Right_arg_value_type const right_scalar;

            //! Constructor.
            explicit Binary_operator_right_scalar(Left_expression const &_left_expression,
                                                  Right_arg_value_type const _right_scalar) :
                    left_expression{_left_expression.cast()}, right_scalar{_right_scalar} {}

            //! Implementation of the \c check_dimensions method, checks the dimensions in the arguments.
            bool check_dimensions_() const { return left_expression.check_dimensions(); }

            //! Implementation of the \c get_size method.
            int get_size_() const { return left_expression.get_size(); }

            //! Implementation of the \c get_dimensions method for the array expression case.
            Dim_array get_dimensions_() const { return left_expression.get_dimensions(); }
        };

        /**
         * Base class for binary operators and functions of two variables helper struct for the instantiation
         * of expression template class.
         * @tparam L left-hand-side argument / first variable type
         * @tparam R right-hand-side argument / second variable type
         * @tparam T return value type
         */
        template<typename L, typename R=L, typename T=typename std::common_type<L, R>::type>
        struct Binary_operator_helper_base {
            using Left_arg_type = L;
            using Right_arg_type = R;
            using Return_type = T;
        };

/**
 * This macro defines an helper utility struct \c Binary_##op_name to be used as template parameter type with
 * the Binary_operator expression template class in order to overload the associated operator \c operator_sign.
 */
#define __EXPRESSION_TEMPLATE_BINARY_OPERATOR_HELPER(op_name, op_sign) \
template<typename L,typename R=L> struct Binary_##op_name : Binary_operator_helper_base<L,R> {\
    using Base = Binary_operator_helper_base<L,R>; \
    static constexpr typename Base::Return_type evaluate(L const l,R const r) {return l op_sign r;} \
};

/**
 * This macro defines an helper utility struct \c Binary_##function to be used as template parameter type with
 * the Binary_operator expression template class in order to overload the corresponding STL mathematic function.
 * \c function must be an STL binary function identifier (without namespace scope resolution operator).
 */
#define __EXPRESSION_TEMPLATE_BINARY_FUNCTION_HELPER(function) \
template<typename L,typename R=L> struct Binary_##function : Binary_operator_helper_base<L,R> {\
    using Base = Binary_operator_helper_base<L,R>; \
    static constexpr typename Base::Return_type evaluate(L const l,R const r) {return std:: function (l,r);} \
};

/**
 * This macro declares and defines all the necessary overloads to extends a STL binary operator or function of two
 * variable to the expression template implementation proposed with Kadath. \c op_name must be the corresponding name
 * identifier used to define the helper utility struct, \c op_sign the operator sign or function STL identifier
 * (without namespace scope resolution operator). \c op_identifier_prefix must respectively be \c operator for a
 * binary arithmetic oeprator or nothing at all for a function of two variables.
 * Note that the \c EXPRESSION_TEMPLATE_BINARY_OPERATOR_OVERLOADS and the
 * \c EXPRESSION_TEMPLATE_BINARY_FUNCTION_OVERLAODS preprocessor macros propose easier to use interfaces, so the
 * following one is an internal one which should not be used.
 */
#define __EXPRESSION_TEMPLATE_BINARY_OPERATOR_OVERLOADS(op_name, op_sign, op_identifier_prefix) \
template<typename L_expr,typename R_expr> using Binary_operator_##op_name##_expr_expr = \
    Binary_operator<\
        Binary_##op_name<typename L_expr::Value_type,typename R_expr::Value_type>,\
        Expression<L_expr>,Expression<R_expr>\
    >; \
template<typename L_expr,typename R_expr> inline Binary_operator_##op_name##_expr_expr<L_expr,R_expr> \
op_identifier_prefix op_sign (Expression<L_expr> const & left_expression,Expression<R_expr> const & right_expression) {\
    return Binary_operator_##op_name##_expr_expr<L_expr,R_expr>{left_expression,right_expression};\
}\
template<typename L,typename R_expr> using Binary_operator_##op_name##_scal_expr = \
    Binary_operator_left_scalar<Binary_##op_name<L,typename R_expr::Value_type>,Expression<R_expr>>;\
template<typename L,typename R_expr> inline Binary_operator_##op_name##_scal_expr<L,R_expr> op_identifier_prefix op_sign \
(L const left_scalar, Expression<R_expr> const & right_expression) {\
    return Binary_operator_##op_name##_scal_expr<L,R_expr>{Scalar_expression<L>{left_scalar},\
                                                           right_expression};\
}\
template<typename L_expr,typename R> using Binary_operator_##op_name##_expr_scal = \
    Binary_operator_right_scalar<Binary_##op_name<typename L_expr::Value_type,R>,Expression<L_expr>>;\
template<typename L_expr,typename R> inline Binary_operator_##op_name##_expr_scal<L_expr,R> op_identifier_prefix op_sign \
(Expression<L_expr> const & left_expression,R const right_scalar) {\
    return Binary_operator_##op_name##_expr_scal<L_expr,R>{left_expression,Scalar_expression<R>{right_scalar}};\
}\
template<typename L,typename R_expr> using Binary_operator_##op_name##_vect_expr = \
    Binary_operator<Binary_##op_name<L,typename R_expr::Value_type>,Vector_expression<L>,Expression<R_expr>>;\
template<typename L,typename R_expr> inline Binary_operator_##op_name##_vect_expr<L,R_expr> op_identifier_prefix op_sign \
(Array<L> const & left_vector, Expression<R_expr> const & right_expression) {\
    return Binary_operator_##op_name##_vect_expr<L,R_expr>{Vector_expression<L>{left_vector},\
                                                           right_expression};\
}\
template<typename L_expr,typename R> using Binary_operator_##op_name##_expr_vect = \
    Binary_operator<Binary_##op_name<typename L_expr::Value_type,R>,Expression<L_expr>,Vector_expression<R>>;\
template<typename L_expr,typename R> inline Binary_operator_##op_name##_expr_vect<L_expr,R> op_identifier_prefix op_sign \
(Expression<L_expr> const & left_expression,Array<R> const & right_vector) {\
    return Binary_operator_##op_name##_expr_vect<L_expr,R>{left_expression,Vector_expression<R>{right_vector}};\
}\
template<typename L,typename R>\
inline Binary_operator<Binary_##op_name <L,R>,Vector_expression<L>,Vector_expression<R>> \
op_identifier_prefix op_sign (Array<L> const & left_vector, Array<R> const & right_vector){\
   return Binary_operator<Binary_##op_name <L,R>,Vector_expression<L>,Vector_expression<R>>{\
                Vector_expression<L>{left_vector},Vector_expression<R>{right_vector}};\
} \
template<typename L,typename R>\
inline Binary_operator_left_scalar<Binary_##op_name<L,R>,Vector_expression<R>> \
op_identifier_prefix op_sign (L const left_scalar,Array<R> const & right_vector) {\
    return Binary_operator_left_scalar<Binary_##op_name<L,R>,Vector_expression<R>>{left_scalar,\
                                                           Vector_expression<R>{right_vector}};\
} \
template<typename L,typename R>\
inline Binary_operator_right_scalar<Binary_##op_name<L,R>,Vector_expression<L>> \
op_identifier_prefix op_sign (Array<L> const & left_vector,R const right_scalar) {\
    return Binary_operator_right_scalar<Binary_##op_name<L,R>,Vector_expression<L>>{\
                Vector_expression<L>{left_vector},right_scalar};\
}

/**
 * Preprocessor macro for the declaration and definition of all the overloads for a single binary
 * arithmetic operator necessary for its use with the expression template optimization.
 * \c operator_name identifier to link the operator to its expression template helper utility class.
 * \c operator_sign sign used to call the operator in an expression (ex: +, -, *, /, etc...).
 */
#define EXPRESSION_TEMPLATE_BINARY_OPERATOR_OVERLOADS(operator_name, operator_sign) \
    __EXPRESSION_TEMPLATE_BINARY_OPERATOR_HELPER(operator_name,operator_sign) \
    __EXPRESSION_TEMPLATE_BINARY_OPERATOR_OVERLOADS(operator_name,operator_sign,operator)

/**
 * Preprocessor macro for the declaration and definition of all the overloads for a single function
 * of two variable necessary for its use within the expression template optimization.
 * \c function identifier of the function to overload.
 */
#define EXPRESSION_TEMPLATE_BINARY_FUNCTION_OVERLOADS(function) \
    __EXPRESSION_TEMPLATE_BINARY_FUNCTION_HELPER(function) \
    __EXPRESSION_TEMPLATE_BINARY_OPERATOR_OVERLOADS(function,function,)

        namespace operators {
            // Overloads for common operators and two variables functions.
            EXPRESSION_TEMPLATE_BINARY_OPERATOR_OVERLOADS(add, +);
            EXPRESSION_TEMPLATE_BINARY_OPERATOR_OVERLOADS(sub, -);
            EXPRESSION_TEMPLATE_BINARY_OPERATOR_OVERLOADS(mul, *);
            EXPRESSION_TEMPLATE_BINARY_OPERATOR_OVERLOADS(div, /);
            EXPRESSION_TEMPLATE_BINARY_FUNCTION_OVERLOADS(pow);
            EXPRESSION_TEMPLATE_BINARY_FUNCTION_OVERLOADS(atan2);
        } // namespace operators
    } // namespace internal
} //namespace Kadath

#endif //__EXPRESSION_BINARY_OPERATORS_HPP_
