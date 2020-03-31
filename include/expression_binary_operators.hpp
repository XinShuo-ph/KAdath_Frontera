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
            Left_expression_type const left_expression;
            //! Right expression sub-expression tree.
            Right_expression_type const right_expression;

            //! Constructor.
            explicit Binary_operator(Left_expression const &_left_expression, Right_expression const &_right_expression) :
                    left_expression{_left_expression.cast()}, right_expression{_right_expression.cast()} {}

            //! Call to the actual evaluation of the operator.
            Value_type evaluate_(int i) const {
                return Operator::evaluate(left_expression.evaluate(i),right_expression.evaluate(i));
            }

            //! Implementation of the \c check_dimensions method, checks the dimensions in the arguments.
            bool check_dimensions_() const {
                return left_expression.check_dimensions() && right_expression.check_dimensions() &&
                       left_expression.get_dimensions() == right_expression.get_dimensions();
            }

            //! Implementation of the \c get_size method.
            int get_size_() const { return left_expression.get_size(); }

            //! Implementation of the \c get_dimensions method for the array expression case.
            Dim_array get_dimensions_() const { return left_expression.get_dimensions(); }
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

            static_assert(std::is_convertible<Right_arg_value_type,Value_type>::value &&
                            std::is_convertible<typename Right_expression::Value_type,Value_type>::value);

            //! Left scalar value.
            Left_arg_value_type const left_scalar;
            //! Right expression sub-expression tree.
            Right_expression_type const right_expression;

            //! Constructor.
            explicit Binary_operator_left_scalar(Left_arg_value_type const _left_scalar,
                                                 Right_expression const &_right_expression) :
                    left_scalar{_left_scalar}, right_expression{_right_expression.cast()} {}

            //! Call to the actual evaluation of the operator.
            Value_type evaluate_(int i) const {return Operator::evaluate(left_scalar,right_expression.evaluate(i));}

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

            static_assert(std::is_convertible<Right_arg_value_type,Value_type>::value &&
                          std::is_convertible<typename Left_expression::Value_type,Value_type>::value);

            //! Right expression sub-expression tree.
            Left_expression_type const left_expression;
            //! Right scalar value.
            Right_arg_value_type const right_scalar;

            //! Constructor.
            explicit Binary_operator_right_scalar(Left_expression const &_left_expression,
                                                  Right_arg_value_type const _right_scalar) :
                    left_expression{_left_expression.cast()}, right_scalar{_right_scalar} {}


            //! Call to the actual evaluation of the operator.
            Value_type evaluate_(int i) const {return Operator::evaluate(left_expression.evaluate(i),right_scalar);}

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
 * This macro declares and defines all the necessary overloads to extends an arithmetic binary operator.
 * \c op_name defines an identifier to link the overloads to the actual operator implementation, any identifier can
 * be use with the constraint that each can only be used once (or different struct with the same name will be declared)
 * \c op_symbol the operator symbol.
 */
#define EXPRESSION_TEMPLATE_BINARY_OPERATOR_OVERLOADS(op_name,op_symbol) \
template<typename L,typename R=L> struct Binary_##op_name : Binary_operator_helper_base<L,R> {\
    using Base = Binary_operator_helper_base<L,R>; \
    static constexpr typename Base::Return_type evaluate(L const l,R const r) {return l op_symbol r;} \
}; \
template<typename L_expr,typename R_expr> using Binary_operator_##op_name##_expr_expr = \
    Binary_operator<\
        Binary_##op_name<typename L_expr::Value_type,typename R_expr::Value_type>,\
        Expression<L_expr>,Expression<R_expr>\
    >; \
template<typename L_expr,typename R_expr> inline Binary_operator_##op_name##_expr_expr<L_expr,R_expr> \
operator op_symbol (Expression<L_expr> const & left_expression,Expression<R_expr> const & right_expression) {\
    return Binary_operator_##op_name##_expr_expr<L_expr,R_expr>{left_expression,right_expression};\
}\
template<typename L,typename R_expr> using Binary_operator_##op_name##_scal_expr = \
    Binary_operator_left_scalar<Binary_##op_name<L,typename R_expr::Value_type>,Expression<R_expr>>;\
template<typename L,typename R_expr> inline typename \
std::enable_if<std::is_arithmetic<L>::value,Binary_operator_##op_name##_scal_expr<L,R_expr>>::type \
operator op_symbol (L const left_scalar, Expression<R_expr> const & right_expression) {\
    return Binary_operator_##op_name##_scal_expr<L,R_expr>{Scalar_expression<L>{left_scalar},\
                                                           right_expression};\
}\
template<typename L_expr,typename R> using Binary_operator_##op_name##_expr_scal = \
    Binary_operator_right_scalar<Binary_##op_name<typename L_expr::Value_type,R>,Expression<L_expr>>;\
template<typename L_expr,typename R> inline typename \
std::enable_if<std::is_arithmetic<R>::value,Binary_operator_##op_name##_expr_scal<L_expr,R>>::type \
operator op_symbol (Expression<L_expr> const & left_expression,R const right_scalar) {\
    return Binary_operator_##op_name##_expr_scal<L_expr,R>{left_expression,Scalar_expression<R>{right_scalar}};\
}\
template<typename L,typename R_expr> using Binary_operator_##op_name##_vect_expr = \
    Binary_operator<Binary_##op_name<L,typename R_expr::Value_type>,Vector_expression<L>,Expression<R_expr>>;\
template<typename L,typename R_expr> inline Binary_operator_##op_name##_vect_expr<L,R_expr> operator op_symbol \
(Array<L> const & left_vector, Expression<R_expr> const & right_expression) {\
    return Binary_operator_##op_name##_vect_expr<L,R_expr>{Vector_expression<L>{left_vector},\
                                                           right_expression};\
}\
template<typename L_expr,typename R> using Binary_operator_##op_name##_expr_vect = \
    Binary_operator<Binary_##op_name<typename L_expr::Value_type,R>,Expression<L_expr>,Vector_expression<R>>;\
template<typename L_expr,typename R> inline Binary_operator_##op_name##_expr_vect<L_expr,R> operator op_symbol \
(Expression<L_expr> const & left_expression,Array<R> const & right_vector) {\
    return Binary_operator_##op_name##_expr_vect<L_expr,R>{left_expression,Vector_expression<R>{right_vector}};\
}\
template<typename L,typename R>\
inline Binary_operator<Binary_##op_name <L,R>,Vector_expression<L>,Vector_expression<R>> \
operator op_symbol (Array<L> const & left_vector, Array<R> const & right_vector){\
   return Binary_operator<Binary_##op_name <L,R>,Vector_expression<L>,Vector_expression<R>>{\
                Vector_expression<L>{left_vector},Vector_expression<R>{right_vector}};\
} \
template<typename L,typename R> inline typename \
std::enable_if<std::is_arithmetic<L>::value,\
               Binary_operator_left_scalar<Binary_##op_name<L,R>,Vector_expression<R>>>::type \
operator op_symbol (L const left_scalar,Array<R> const & right_vector) {\
    return Binary_operator_left_scalar<Binary_##op_name<L,R>,Vector_expression<R>>{left_scalar,\
                                                           Vector_expression<R>{right_vector}};\
} \
template<typename L,typename R> inline typename \
std::enable_if<std::is_arithmetic<R>::value, \
               Binary_operator_right_scalar<Binary_##op_name<L,R>,Vector_expression<L>>>::type \
operator op_symbol (Array<L> const & left_vector,R const right_scalar) {\
    return Binary_operator_right_scalar<Binary_##op_name<L,R>,Vector_expression<L>>{\
                Vector_expression<L>{left_vector},right_scalar};\
}

/**
 * This macro declares and defines all the necessary overloads to extends an STL function of two variables.
 * \c function is the function identifier without its eventual scope resolution operator
 */
#define EXPRESSION_TEMPLATE_BINARY_FUNCTION_OVERLOADS(function) \
template<typename L,typename R=L> struct Binary_##function : Binary_operator_helper_base<L,R> {\
    using Base = Binary_operator_helper_base<L,R>; \
    static constexpr typename Base::Return_type evaluate(L const l,R const r) {return std:: function (l,r);} \
}; \
template<typename L_expr,typename R_expr> using Binary_operator_##function##_expr_expr = \
    Binary_operator<\
        Binary_##function<typename L_expr::Value_type,typename R_expr::Value_type>,\
        Expression<L_expr>,Expression<R_expr>\
    >; \
template<typename L_expr,typename R_expr> inline Binary_operator_##function##_expr_expr<L_expr,R_expr> \
function (Expression<L_expr> const & left_expression,Expression<R_expr> const & right_expression) {\
    return Binary_operator_##function##_expr_expr<L_expr,R_expr>{left_expression,right_expression};\
}\
template<typename L,typename R_expr> using Binary_operator_##function##_scal_expr = \
    Binary_operator_left_scalar<Binary_##function<L,typename R_expr::Value_type>,Expression<R_expr>>;\
template<typename L,typename R_expr> inline typename \
std::enable_if<std::is_arithmetic<L>::value,Binary_operator_##function##_scal_expr<L,R_expr>>::type \
function (L const left_scalar, Expression<R_expr> const & right_expression) {\
    return Binary_operator_##function##_scal_expr<L,R_expr>{Scalar_expression<L>{left_scalar},\
                                                           right_expression};\
}\
template<typename L_expr,typename R> using Binary_operator_##function##_expr_scal = \
    Binary_operator_right_scalar<Binary_##function<typename L_expr::Value_type,R>,Expression<L_expr>>;\
template<typename L_expr,typename R> inline typename \
std::enable_if<std::is_arithmetic<R>::value,Binary_operator_##function##_expr_scal<L_expr,R>>::type \
function (Expression<L_expr> const & left_expression,R const right_scalar) {\
    return Binary_operator_##function##_expr_scal<L_expr,R>{left_expression,Scalar_expression<R>{right_scalar}};\
}\
template<typename L,typename R_expr> using Binary_operator_##function##_vect_expr = \
    Binary_operator<Binary_##function<L,typename R_expr::Value_type>,Vector_expression<L>,Expression<R_expr>>;\
template<typename L,typename R_expr> inline Binary_operator_##function##_vect_expr<L,R_expr> function \
(Array<L> const & left_vector, Expression<R_expr> const & right_expression) {\
    return Binary_operator_##function##_vect_expr<L,R_expr>{Vector_expression<L>{left_vector},\
                                                           right_expression};\
}\
template<typename L_expr,typename R> using Binary_operator_##function##_expr_vect = \
    Binary_operator<Binary_##function<typename L_expr::Value_type,R>,Expression<L_expr>,Vector_expression<R>>;\
template<typename L_expr,typename R> inline Binary_operator_##function##_expr_vect<L_expr,R> function \
(Expression<L_expr> const & left_expression,Array<R> const & right_vector) {\
    return Binary_operator_##function##_expr_vect<L_expr,R>{left_expression,Vector_expression<R>{right_vector}};\
}\
template<typename L,typename R>\
inline Binary_operator<Binary_##function <L,R>,Vector_expression<L>,Vector_expression<R>> \
function (Array<L> const & left_vector, Array<R> const & right_vector){\
   return Binary_operator<Binary_##function <L,R>,Vector_expression<L>,Vector_expression<R>>{\
                Vector_expression<L>{left_vector},Vector_expression<R>{right_vector}};\
} \
template<typename L,typename R> inline typename \
std::enable_if<std::is_arithmetic<L>::value,\
               Binary_operator_left_scalar<Binary_##function<L,R>,Vector_expression<R>>>::type \
function (L const left_scalar,Array<R> const & right_vector) {\
    return Binary_operator_left_scalar<Binary_##function<L,R>,Vector_expression<R>>{left_scalar,\
                                                           Vector_expression<R>{right_vector}};\
} \
template<typename L,typename R> inline typename \
std::enable_if<std::is_arithmetic<R>::value,\
               Binary_operator_right_scalar<Binary_##function<L,R>,Vector_expression<L>>>::type \
function (Array<L> const & left_vector,R const right_scalar) {\
    return Binary_operator_right_scalar<Binary_##function<L,R>,Vector_expression<L>>{\
                Vector_expression<L>{left_vector},right_scalar};\
}


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
