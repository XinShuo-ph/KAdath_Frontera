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
#include "array.hpp"
#include "index.hpp"

using namespace Kadath;

int main(int argc,char *argv[])
{
    constexpr int const ndim{5};
    std::array<int,ndim> d={10,9,8,11,13};
    Dim_array dimensions{ndim};
    dimensions.set(0) = d[0];
    dimensions.set(1) = d[1];
    dimensions.set(2) = d[2];
    dimensions.set(3) = d[3];
    dimensions.set(4) = d[4];
    int const nbr {d[0]*d[1]*d[2]*d[3]*d[4]};
    Array<int> test_array{dimensions};
    for(int i=0;i<nbr;i++) test_array.set_data()[i] = i;
    Index index{dimensions};
    Array_iterator array_index{dimensions};
    {
        bool inc_ok{true};
        bool inc_i{true};
        bool inc_ai{true};
        do {
            int const v1{test_array(index)};
            int const v2{test_array(array_index)};
            EXPECT_EQ(v1, v2);
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
                EXPECT_EQ(v1, v2);
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
    for(int var=0;var<ndim;var++)
        for(int step=2;step<d[var];step++) {
            index.set_start();
            array_index.set_start();
            bool incstep_ok{true};
            do {
                int const v1{test_array(index)};
                int const v2{test_array(array_index)};
                EXPECT_EQ(v1, v2);
                incstep_ok = v1 == v2;
                inc_i = index.inc(step, var);
                inc_ai = array_index.inc(step, var);
            } while (incstep_ok && inc_i && inc_ai);
            assert(incstep_ok);
        }
    return 0;
}
