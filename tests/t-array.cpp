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

/**
 * \file t-array.cpp
 * This sets of unit tests checks the functions used to switch multi-dimensional
 * array ordering .
 */

#include <random>
#include <string>
#include "array.hpp"

void save_array(Kadath::Array<int> const & array,std::string const &file_name,
                Kadath::Array_ordering order = Kadath::last_index);
Kadath::Array<int> read_array(std::string const & file_name,
                              Kadath::Array_ordering order = Kadath::last_index);
bool operator==(Kadath::Array<int> const &l,Kadath::Array<int> const &r);

class Array_tester {
public:
    static constexpr int one_dim_array_max_size {1000};
    static constexpr int two_dim_array_max_size {100};
    static constexpr int multi_dim_array_max_size {10};
    static constexpr int int_max {std::numeric_limits<int>::max()};
private:
    std::uniform_int_distribution<> distribution;
    std::mt19937 generator;

public:
    Array_tester(int a=0,int b=int_max) : distribution{a,b}, generator{std::random_device{}()} {}

    bool test_1_1();
    bool test_1_x(unsigned nb_check = 1);
    bool test_2_x(unsigned nb_check = 1);
    bool test_x_x(unsigned nb_check=1,int ndim_max = 5) ;
    bool test_file_1_x(unsigned nb_check = 1);
    bool test_file_2_x(unsigned nb_check = 1);
    bool test_file_x_x(unsigned nb_check = 1,int ndim_max = 5);

    int random() {return distribution(generator);}
    int random(int a, int b=int_max) {return distribution(generator, std::uniform_int_distribution<>::param_type{a, b});}
    void random_fill(Kadath::Array<int> & array,int min_value = 0,
                     int max_value = int_max)
    {for(int k=0;k<array.get_nbr();k++) array.set_data()[k] = this->random(min_value,max_value);};
    void sequence_fill(Kadath::Array<int> & array,int start = 0,int step=1)
    {for(int k=0;k<array.get_nbr();k++) array.set_data()[k] = start + k*step;}
};

int main(int argc,char * argv[]) {
    Array_tester array_tester{};
    std::cout << "========== t-array unit-tests set ==========\n\n";

    std::cout << "Tests for reordering methods : " << std::endl;
    std::cout << "    - trivial zero-dimensional array case (one element array)...";
    bool const test_1_1 {array_tester.test_1_1()};
    std::cout << (test_1_1 ? " success" : " failure") << std::endl;
    assert(test_1_1);

    std::cout << "    - trivial 1-dimensional array case...";
    bool const test_1_x {array_tester.test_1_x(100)};
    std::cout << (test_1_x ? " success" : " failure") << std::endl;
    assert(test_1_x);

    std::cout << "    - 2-dimensional array case...";
    bool const test_2_x {array_tester.test_2_x(100)};
    std::cout << (test_2_x ? " success" : " failure") << std::endl;
    assert(test_2_x);

    std::cout << "    - higher-dimensional array case...";
    bool const test_x_x {array_tester.test_x_x(10)};
    std::cout << (test_x_x ? " success" : " failure") << std::endl;
    assert(test_x_x);

    std::cout << std::endl <<
        "Tests for file I/O methods mixed with reordering :" << std::endl;
    std::cout << "    - trivial 1-dimensional array case...";
    bool const test_file_1_x {array_tester.test_file_1_x(10)};
    std::cout << (test_file_1_x ? " success" : " failure") << std::endl;
    assert(test_file_1_x);

    std::cout << "    - 2-dimensional array case...";
    bool const test_file_2_x {array_tester.test_file_2_x(10)};
    std::cout << (test_file_2_x ? " success" : " failure") << std::endl;
    assert(test_file_2_x);

    std::cout << "    - higher-dimensional array case";
    bool const test_file_x_x {array_tester.test_file_x_x(5)};
    std::cout << (test_file_x_x ? " success" : " failure") << std::endl;
    assert(test_file_x_x);

    std::cout << "\n\n=============================================";
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

bool Array_tester::test_file_1_x(unsigned int nb_check) {
    bool success {true};
    for(int n=0;n<nb_check && success;n++)
    {
        std::string const fname_li{"t-array_1_x_li.dat"}, fname_fi{"t-array_1_x_fi.dat"};
        Kadath::Array<int> array{random(1,one_dim_array_max_size)};
        random_fill(array);
        save_array(array,fname_li);
        save_array(array,fname_fi,Kadath::first_index);
        Kadath::Array<int> array_li{read_array(fname_li)},
                            array_fi{read_array(fname_fi,Kadath::first_index)};
        bool const li_ok{array_li == array};
        bool const fi_ok{array_fi == array};
        success = (li_ok && fi_ok);
    }
    return success;
}

bool Array_tester::test_file_2_x(unsigned int nb_check) {
    bool success {true};
    for(int n=0;n<nb_check && success;n++)
    {
        std::string const fname_li{"t-array_2_x_li.dat"}, fname_fi{"t-array_2_x_fi.dat"};
        int const M{random(1,two_dim_array_max_size)},
                  N{random(1,two_dim_array_max_size)};
        Kadath::Array<int> array{M,N};
        random_fill(array);
        save_array(array,fname_li);
        save_array(array,fname_fi,Kadath::first_index);
        Kadath::Array<int> array_li{read_array(fname_li)},
                array_fi{read_array(fname_fi,Kadath::first_index)};
        bool const li_ok{array_li == array};
        bool const fi_ok{array_fi == array};
        success = (li_ok && fi_ok);
    }
    return success;
}

bool Array_tester::test_file_x_x(unsigned int nb_check,int ndim_max) {
    bool success {true};
    for(int n=0;n<nb_check && success;n++)
    {
        std::string const fname_li{"t-array_x_x_li.dat"}, fname_fi{"t-array_x_x_fi.dat"};
        int const ndim {random(3,ndim_max)};
        Kadath::Dim_array dimensions {ndim};
        for(int i=0;i<ndim;i++) dimensions.set(i) = random(1,multi_dim_array_max_size);
        Kadath::Array<int> array{dimensions};
        random_fill(array);
        save_array(array,fname_li);
        save_array(array,fname_fi,Kadath::first_index);
        Kadath::Array<int> array_li{read_array(fname_li)},
                array_fi{read_array(fname_fi,Kadath::first_index)};
        bool const li_ok{array_li == array};
        bool const fi_ok{array_fi == array};
        success = (li_ok && fi_ok);
    }
    return success;
}

void save_array(Kadath::Array<int> const &array,std::string const &file_name,
                Kadath::Array_ordering order) {
    FILE *file = fopen(file_name.c_str(), "w");
    array.save(file,order);
    fclose(file);
}

Kadath::Array<int> read_array(std::string const & file_name,
                              Kadath::Array_ordering order) {
    FILE *file = fopen(file_name.c_str(),"r");
    Kadath::Array<int> array{file,order};
    fclose(file);
    return std::move(array);
}

bool operator==(Kadath::Array<int> const &l,Kadath::Array<int> const &r) {
    bool success{l.get_dimensions() == r.get_dimensions()};
    for(int i=0;i<l.get_nbr() && success;i++)
        success = (l.get_data()[i] == r.get_data()[i]);
    return success;
}