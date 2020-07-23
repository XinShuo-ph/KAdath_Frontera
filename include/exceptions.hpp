//
// Created by sauliac on 29/06/2020.
//

#ifndef __EXCEPTIONS_HPP_
#define __EXCEPTIONS_HPP_

#include "headcpp.hpp"
#include <exception>

namespace Kadath {


class Unknown_base_error : public std::runtime_error {
public:
    int base;
    std::string const where;
    Unknown_base_error(int b,std::string const & w) : std::runtime_error{"Unknown base"},
        base{b}, where{w} {}

    const char * what() const noexcept override {
        std::string explanation {"Unknown base in"};
        explanation += where + " - base code = ";
        explanation += std::to_string(base);
        return explanation.c_str();
    }
};


}

#endif //__EXCEPTIONS_HPP_
