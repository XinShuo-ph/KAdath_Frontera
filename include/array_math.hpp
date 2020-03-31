/*
    Copyright 2017 Philippe Grandclement

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

#include "expression_unary_operators.hpp"
#include "expression_binary_operators.hpp"

namespace Kadath {

    // put all expression template operators and functions overloads declarations in the Kadath namespace
    using namespace internal::operators;


    template<typename T> template<typename E> Array<T>::Array(internal::Expression<E> const & expression) :
            dimensions{expression.get_dimensions()}, nbr{expression.get_size()}, data{new T[nbr]}
    {
        static_assert(std::is_convertible<typename E::Value_type,T>::value,"Expression value type must be convertible"
                                                                           " to the target underlying type.");
        assert(expression.check_dimensions() && dimensions == expression.get_dimensions());
        for(int i{0};i<nbr;i++) data[i] = expression.evaluate(i);
    }

    template<typename T> template<typename E> Array<T>& Array<T>::operator=(internal::Expression<E> const & expression)
    {
        static_assert(std::is_convertible<typename E::Value_type,T>::value,"Expression value type must be convertible"
                                                                           " to the target underlying type.");
        assert(expression.check_dimensions() && dimensions == expression.get_dimensions());
        for(int i{0};i<nbr;i++) data[i] = expression.evaluate(i);
    }

    template<typename T> template<typename Expr> Array<T>& Array<T>::operator+=(internal::Expression<Expr> const & expression)
    {
        static_assert(std::is_convertible<typename Expr::Value_type,T>::value,"Expression value type must be convertible"
                                                                              " to the target underlying type.");
        assert(expression.check_dimensions() && dimensions == expression.get_dimensions());
        for(int i{0};i<nbr;i++) data[i] += expression.evaluate(i);
        return *this;
    }

    template<typename T> template<typename Expr> Array<T>& Array<T>::operator-=(internal::Expression<Expr> const & expression)
    {
        static_assert(std::is_convertible<typename Expr::Value_type,T>::value,"Expression value type must be convertible"
                                                                              " to the target underlying type.");
        assert(expression.check_dimensions() && dimensions == expression.get_dimensions());
        for(int i{0};i<nbr;i++) data[i] -= expression.evaluate(i);
        return *this;
    }

    template<typename T> template<typename Expr> Array<T>& Array<T>::operator*=(internal::Expression<Expr> const & expression)
    {
        static_assert(std::is_convertible<typename Expr::Value_type,T>::value,"Expression value type must be convertible"
                                                                              " to the target underlying type.");
        assert(expression.check_dimensions() && dimensions == expression.get_dimensions());
        for(int i{0};i<nbr;i++) data[i] *= expression.evaluate(i);
        return *this;
    }

    template<typename T> template<typename Expr> Array<T>& Array<T>::operator/=(internal::Expression<Expr> const & expression)
    {
        static_assert(std::is_convertible<typename Expr::Value_type,T>::value,"Expression value type must be convertible"
                                                                              " to the target underlying type.");
        assert(expression.check_dimensions() && dimensions == expression.get_dimensions());
        for(int i{0};i<nbr;i++) data[i] /= expression.evaluate(i);
        return *this;
    }



template <typename T> inline T scal (const Array<T>& a, const Array<T>& b) {
	T res = a.data[0] * b.data[0] ;
	for (int i=1 ; i<a.get_size(0) ; i++)
		res += a.data[i] * b.data[i] ;
	return res ;
}


template <typename T> inline T diffmax (const Array<T>& a, const Array<T>& b) {
	assert (a.nbr==b.nbr) ;
	T res = 0 ;
	T diff ;
	for (int i=0 ; i<a.nbr ; i++) {
		diff = fabs(a.data[i]-b.data[i]) ;
		if (diff> res) 
			res = diff ;
	}
	return res ;
}

template <typename T> inline  T max (const Array<T>& so)  {
	T res = so.data[0] ;
	for (int i=1 ; i<so.nbr ; i++)
		if (so.data[i]>res)
			res = so.data[i] ;
	return res;
}

template<typename Expr> inline typename Expr::Value_type max(internal::Expression<Expr> const &expression)
{
    typename Expr::Value_type res {expression.evaluate(0)};
    int const size {expression.get_size()};
    for(int i{1};i<size;i++) res = std::max(expression.evaluate(i),res);
    return res;
}

template <typename T> inline  T min (const Array<T>& so)  {
	T res = so.data[0] ;
	for (int i=1 ; i<so.nbr ; i++)
		if (so.data[i]<res)
			res = so.data[i] ;
	return res;
}

template<typename Expr> inline typename Expr::Value_type min(internal::Expression<Expr> const &expression)
{
    typename Expr::Value_type res {expression.evaluate(0)};
    int const size {expression.get_size()};
    for(int i{1};i<size;i++) res = std::min(expression.evaluate(i),res);
    return res;
}

template <typename T> inline  T sum (const Array<T>& so)  {
	T res = 0 ;
	for (int i=0 ; i<so.nbr ; i++)
		res += so.data[i] ;
	return res;
}

template<typename Expr> inline typename Expr::Value_type sum(internal::Expression<Expr> const &expression)
{
    typename Expr::Value_type res {0};
    int const size {expression.get_size()};
    for(int i{0};i<size;i++) res += expression.evaluate(i);
    return res;
}
}
