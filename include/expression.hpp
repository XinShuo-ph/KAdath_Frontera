//
// Created by sauliac on 27/03/2020.
//

#include <type_traits>
#include <array.hpp>
#include <cassert>
#include <exception>

#ifndef __EXPRESSION_HPP_
#define __EXPRESSION_HPP_

namespace Kadath {

    namespace internal {

        /**
         * Dummy void-wrapping type used for template parameter type or type alias when a type is unknown.
         */
        struct Undefined_type {
            using type = void;
            static constexpr bool value {false};
            Undefined_type() = delete;
        };


        /**
         * Base class for all expressions.
         * @tparam E actual type of the expression.
         * @tparam T underlying value type.
         */
        template<typename E,typename T> struct Expression {
            //! Type alias for the actual expression type.
            using Expression_type = E;
            //! Type alias for the underlying value type.
            using Value_type = T;

            //! Casts the expression to its true type.
            Expression_type const & cast() const {return static_cast<Expression_type const &>(*this);}
            /**
             * Evaluate the expression for the \c i-th component.
             * @param i index of the component to evaluate.
             * @return the value of the result of the expression at index \c i.
             */
            Value_type evaluate(int i) const {return this->cast().evaluate_(i);}
            //! Checks if the dimensions of the array leaves of the expression tree are all equals.
            bool check_dimensions() const {return this->cast().check_dimensions_();}
            //! Returns the size of the first array evaluated when building the expression tree.
            int get_size() const {return this->cast().get_size_();}
            //! Returns the \c dimensions member of the first array evaluated while building the expression tree.
            Dim_array get_dimensions() const {return this->cast().get_dimensions_();}
        };

        /**
         * Expression type for scalar constants.
         * @tparam T type of the scalar constant.
         */
        template<typename T> struct Scalar_expression : public Expression<Scalar_expression<T>,T> {
            //! Type alias for the type of the constant.
            using Value_type = T;
            //! Alias toward the actual expression type.
            using Expression_type = Scalar_expression<T>;

            //! Value of the scalar constant.
            Value_type const value;

            //! Constructor.
            explicit Scalar_expression(Value_type const _value) : value{_value} {}

            //! Implementation of the \c evaluate method. Returns the constant's value whatever the index value is.
            Value_type evaluate_(int) const {return value;}
            //! Implementation of the \c check_dimensions method, return true since scalar do not raise dimensional matching issues.
            static constexpr bool check_dimensions_() {return true;}
            //! Implementation of the \c get_size method, this should not be called.
            static int get_size_() {
                std::cerr << "Error: cannot call Scalar_constant::get_size_()" << std::endl;
                throw std::runtime_error{"cannot call Scalar_constant::get_size_()"};
            }
            static Dim_array get_dimensions_() {
                std::cerr << "Error: cannot call Scalar_constant::get_dimensions_()" << std::endl;
                throw std::runtime_error{"cannot call Scalar_constant::get_dimensions_()"};
            }
        };

        /**
         * Expression type for the \c Array class.
         * @tparam T underlying element types.
         */
        template<typename T> struct Vector_expression : public Expression<Vector_expression<T>,T> {
            //! Type alias for the array elements value type.
            using Value_type = T;
            //! Type alias for the actual expression type.
            using Expression_type = Vector_expression<T>;

            //! Reference toward the array instance used in the expression.
            Array<Value_type> const & vector;

            //! Constructor.
            explicit Vector_expression(Array<Value_type> const & _vector) : vector{_vector} {}
            //! Implementation of the \c evaluate method. Returns the \c i-th component of \c vector.
            Value_type evaluate_(int i) const {return vector.get_data()[i];}
            //! Implementation of the \c check_dimensions method, return true since a single vector has no problem.
            static constexpr bool check_dimensions_() {return true;}
            //! Implementation of the \c get_size method.
            int get_size_() const {return vector.get_nbr();}
            //! Implementation of the \c get_dimensions method for the array expression case.
            Dim_array get_dimensions_() const {return vector.get_dimensions();}
        };

    }

}

#endif //__EXPRESSION_HPP_
