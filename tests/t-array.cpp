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

#include <random>
#include "array.hpp"

class Array_tester {
public:
    static constexpr int one_dim_array_max_size {1000};
    static constexpr int two_dim_array_max_size {100};
    static constexpr int multi_dim_array_max_size {10};
private:
    std::uniform_int_distribution<> distribution;
    std::mt19937 generator;

public:
    Array_tester(int a=0,int b=std::numeric_limits<int>::max()) : distribution{a,b}, generator{std::random_device{}()} {}

    bool test_1_1();
    bool test_1_x(unsigned nb_check = 1);
    bool test_2_x(unsigned nb_check = 1);
    bool test_x_x(unsigned nb_check=1,int ndim_max = 5) ;

    int random() {return distribution(generator);}
    int random(int a, int b) {return distribution(generator, std::uniform_int_distribution<>::param_type{a, b});}

};

int main(int argc,char * argv[]) {
    Array_tester array_tester{};

    bool const test_1_1 {array_tester.test_1_1()};
    assert(test_1_1);

    bool const test_1_x {array_tester.test_1_x(100)};
    assert(test_1_x);

    bool const test_2_x {array_tester.test_2_x(100)};
    assert(test_2_x);

    bool const test_x_x {array_tester.test_x_x(10)};
    assert(test_x_x);

    return 0;
}

bool Array_tester::test_1_1() {
    Kadath::Array<int> t{1};
    bool const cmo2rmo {t.to_last_dim_major_index(0) == 0};
    bool const rmo2cmo {t.to_first_dim_major_index(0) == 0};
    return cmo2rmo && rmo2cmo;
}

bool Array_tester::test_1_x(unsigned nb_check) {
    bool success {true};
    for(int n=0;n<nb_check && success;n++) {
        Kadath::Array<int> t{random(1,one_dim_array_max_size)};
        for(int i=0;i<t.get_nbr() && success;i++) {
            bool const cmo2rmo {t.to_last_dim_major_index(i) == i};
            bool const rmo2cmo {t.to_first_dim_major_index(i) == i};
            success = cmo2rmo && rmo2cmo;
        }
    }
    return success;
}

bool Array_tester::test_2_x(unsigned int nb_check) {
    bool success {true};
    auto rnd2d = [this]() -> int {return random(1,two_dim_array_max_size);};
    for(int n=0;n<nb_check && success;n++) {
        int const M{rnd2d()}, N{rnd2d()};
        Kadath::Array<int> t{M,N};
        for(int k=0;k<t.get_nbr();k++) {
            int const i {k % M}, j {k / M};
            int const l {i*N + j};
            success = (t.to_last_dim_major_index(k) == l);
            int const ii {k % N}, jj {k / N};
            int const ll {i + j*M};
            success = (success && t.to_first_dim_major_index(k) == ll);
        }
    }
    return success;
}

bool Array_tester::test_x_x(unsigned int nb_check, int ndim_max) {
    bool success {true};
    for(int n=0;n<nb_check && success;n++) {
        int const ndim {random(3,ndim_max)};
        Kadath::Dim_array dimensions {ndim};
        for(int i=0;i<ndim;i++) dimensions.set(i) = random(1,multi_dim_array_max_size);
        Kadath::Array<int> t{dimensions};
        for(int i=0;i<t.get_nbr() && success;i++) {
            success = (t.to_last_dim_major_index(t.to_first_dim_major_index(i)) == i);
        }
    }
    return success;
}