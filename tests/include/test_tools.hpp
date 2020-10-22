//
// Created by sauliac on 20/10/2020.
//

#include <numeric>
#include <random>
#include <type_traits>

#ifndef __KADATH_TESTS_TEST_TOOLS_HPP_
#define __KADATH_TESTS_TEST_TOOLS_HPP_

namespace tests {
    template<typename T> class Random_generator {
    public:
        static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value,
                "Error, T must either be an integral or a floating-point built-in type");
        using value_type = T;
        using arg_type = typename std::conditional<sizeof(T) <= sizeof(T*),T,T const &>::type;
        using distribution_type =
                typename std::conditional<
                    std::is_integral<value_type>::value,
                    std::uniform_int_distribution<value_type>,
                    std::uniform_real_distribution<value_type>
                >::type;

        static constexpr value_type min_val = std::numeric_limits<value_type>::min();
        static constexpr value_type max_val = std::numeric_limits<value_type>::max();

    private:
        distribution_type distribution;
        std::mt19937 generator;
        value_type _min_val;
        value_type _max_val;

    public:
        Random_generator(arg_type a=min_val,arg_type b=max_val) : distribution{a,b}, generator{std::random_device{}()},
            _min_val{a}, _max_val{b} {}
        value_type operator()() {return distribution(generator);}
        value_type operator()(arg_type a,arg_type b) {return distribution(generator, std::uniform_int_distribution<>::param_type{a, b});}

    };
}

#endif //__KADATH_TESTS_TEST_TOOLS_HPP_
