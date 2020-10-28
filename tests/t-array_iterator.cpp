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
#include <cassert>
#include "test_tools.hpp"
#include "array.hpp"
#include "index.hpp"

using namespace Kadath;

int main(int argc,char *argv[])
{
    std::cout << "========================= t-array_iterator unit-tests set ===========================\n\n";

    constexpr int max_size {15};
    tests::Random_generator<int> dist{1,max_size};
    constexpr int const ndim{5};
    std::array<int,ndim> d;
    for(auto & d_i : d) d_i = dist();
    Dim_array dimensions{ndim};
    for(int i{0};i<ndim;i++) dimensions.set(i) = d[i];
    int  nbr {1};
    for(auto d_i : d) nbr *= d_i;
    Array<int> test_array{dimensions};
    for(int i=0;i<nbr;i++) test_array.set_data()[i] = dist(dist.min_val,dist.max_val);

    Index index{dimensions};
    Array_iterator array_index{dimensions};
    {
        bool inc_ok{true};
        bool inc_i{true};
        bool inc_ai{true};
        do {
            int const v1{test_array(index)};
            int const v2{test_array(array_index)};
            assert(v1 == v2);
            inc_ok = v1 == v2;
            inc_i = index.inc();
            inc_ai = array_index.inc();
        } while (inc_ok && inc_ai && inc_i);
        assert(inc_ok);
    }
    {
        bool inc_i{true};
        bool inc_ai{true};
        for (int var = 0; var < ndim; var++) {
            index.set_start();
            array_index.set_start();
            bool inc1_ok{true};
            do {
                int const v1{test_array(index)};
                int const v2{test_array(array_index)};
                assert(v1 == v2);
                inc1_ok = v1 == v2;
                inc_i = index.inc1(var);
                inc_ai = array_index.inc1(var);
            } while (inc1_ok && inc_i && inc_ai);
            assert(inc1_ok);
        }
    }
    {
        bool inc_i{true};
        bool inc_ai{true};
        for (int var = 0; var < ndim; var++)
            for (int step = 2; step < d[var]; step++) {
                index.set_start();
                array_index.set_start();
                bool incstep_ok{true};
                do {
                    int const v1{test_array(index)};
                    int const v2{test_array(array_index)};
                    assert(v1 == v2);
                    incstep_ok = v1 == v2;
                    inc_i = index.inc(step, var);
                    inc_ai = array_index.inc(step, var);
                } while (incstep_ok && inc_i && inc_ai);
                assert(incstep_ok);
            }
    }
    std::cout << "\n\n=====================================================================================";
    return 0;
}
