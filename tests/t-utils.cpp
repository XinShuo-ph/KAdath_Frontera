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
 * Unit tests for utilitary class.
 */

#include <cassert>
#include <type_traits>
#include <thread>
#include "profiled_object.hpp"

using namespace Kadath;

// Large tolerance, because this is no numeric error, but a variation in CPU time measurements of second long
// process, so 1 millisecond is accuracy is good enough. If you encounter trouble with thie tests below,
// 0.1 might still be acceptable.
static constexpr double tol {5.e-3};
template<typename T1,typename T2> inline
typename std::enable_if<std::is_convertible<typename std::common_type<T1,T2>::type,double>::value,bool>::type is_near(T1 const & l,T2 const & r,double _tol = tol)  {
    using T = typename std::common_type<T1,T2>::type;
    T const _l = l, _r = r;
    T const d {std::abs(_l - _r)};
    return d <= _tol;
}

template<typename T> inline bool is_small(T const &x,double _tol = tol) {return is_near(x,0,_tol);}

class AutoProfilerTester : public Profiled_object<AutoProfilerTester> {
public :
    using Base = Profiled_object<AutoProfilerTester>;
    using Base::Base;

    void test_measurement();
    void test_time_conversions();
};

void AutoProfilerTester::test_measurement()
{
    constexpr unsigned short ncheck {10u};
    constexpr unsigned long small_sleep_duration {100}; // milliseconds
    Hash_key long_sleep_key = start_chrono("one_long_sleep"),short_sleeps_key{};
    for(int i=0;i<ncheck;++i)
    {
        Hash_key new_short_sleeps_key = start_chrono("several_short_sleeps");
        if(i==0)
        {
            short_sleeps_key = new_short_sleeps_key;
        }
        else
        {
            assert(new_short_sleeps_key == short_sleeps_key);
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(small_sleep_duration));
        stop_chrono(short_sleeps_key);
    }
    stop_chrono(long_sleep_key);

    this->finalize_profiling();

    auto long_sleep_it = get_statistic_map().find(long_sleep_key);
    assert(long_sleep_it != get_statistic_map().end());
    Statistics const & long_sleep_duration_measure {long_sleep_it->second};
    assert(long_sleep_duration_measure.n_samples==1);
    assert(long_sleep_duration_measure.total_duration >= (ncheck*small_sleep_duration/1000.));
    assert(long_sleep_duration_measure.total_duration == long_sleep_duration_measure.average_duration);
    assert(long_sleep_duration_measure.std_deviation == 0.);
    auto short_sleeps_it = get_statistic_map().find(short_sleeps_key);
    assert(short_sleeps_it != get_statistic_map().end());
    Statistics const & short_sleeps_duration_measure {short_sleeps_it->second};
    assert(short_sleeps_duration_measure.n_samples == ncheck);
    assert(short_sleeps_duration_measure.total_duration >= small_sleep_duration/1000.);
    assert(is_near(short_sleeps_duration_measure.average_duration,small_sleep_duration/1000.));
    assert(is_near(short_sleeps_duration_measure.std_deviation,0.));
}


void AutoProfilerTester::test_time_conversions()
{
    std::chrono::hours hour{1};
    std::chrono::seconds zero{0};
    auto hour_to_hours = to_hours<unsigned>(hour);
    auto hour_to_minutes = to_minutes<unsigned>(hour);
    auto hour_to_seconds = to_seconds<unsigned>(hour);
    auto hour_to_milliseconds = to_milliseconds<unsigned>(hour);
    auto hour_to_microseconds = to_microseconds<unsigned>(hour);
    auto hour_to_nanoseconds = to_nanoseconds<unsigned long>(hour);
    assert(hour_to_hours == 1);
    assert(hour_to_minutes == 60);
    assert(hour_to_seconds == 3600);
    assert(hour_to_milliseconds == 3600000);
    assert(hour_to_microseconds == 3600000000);
    assert(hour_to_nanoseconds == 3600000000000);
    auto zero_to_hours = to_hours<unsigned>(zero);
    auto zero_to_minutes = to_minutes<unsigned>(zero);
    auto zero_to_seconds = to_seconds<unsigned>(zero);
    auto zero_to_milliseconds = to_milliseconds<unsigned>(zero);
    auto zero_to_microseconds = to_microseconds<unsigned>(zero);
    auto zero_to_nanoseconds = to_nanoseconds<unsigned>(zero);
    assert(zero_to_hours == 0);
    assert(zero_to_minutes == 0);
    assert(zero_to_seconds == 0);
    assert(zero_to_milliseconds == 0);
    assert(zero_to_microseconds == 0);
    assert(zero_to_nanoseconds == 0);

    std::chrono::nanoseconds one_nano{1};
    std::chrono::seconds f_zero{0};
    auto one_nano_to_hours = to_hours(one_nano);
    auto one_nano_to_minutes = to_minutes(one_nano);
    auto one_nano_to_seconds = to_seconds(one_nano);
    auto one_nano_to_milliseconds = to_milliseconds(one_nano);
    auto one_nano_to_microseconds = to_microseconds(one_nano);
    auto one_nano_to_nanoseconds = to_nanoseconds(one_nano);
    assert(one_nano_to_hours == 0.000000001/3600.);
    assert(one_nano_to_minutes == 0.000000001/60.);
    assert(one_nano_to_seconds == 0.000000001);
    assert(one_nano_to_milliseconds == 0.000001);
    assert(one_nano_to_microseconds == 0.001);
    assert(one_nano_to_nanoseconds == 1.);
    auto f_zero_to_hours = to_hours(f_zero);
    auto f_zero_to_minutes = to_minutes(f_zero);
    auto f_zero_to_seconds = to_seconds(f_zero);
    auto f_zero_to_milliseconds = to_milliseconds(f_zero);
    auto f_zero_to_microseconds = to_microseconds(f_zero);
    auto f_zero_to_nanoseconds = to_nanoseconds(f_zero);
    assert(f_zero_to_hours == 0.);
    assert(f_zero_to_minutes == 0.);
    assert(f_zero_to_seconds == 0.);
    assert(f_zero_to_milliseconds == 0.);
    assert(f_zero_to_microseconds == 0.);
    assert(f_zero_to_nanoseconds == 0.);

    auto zero_to_hh_mm_ss = to_hh_mm_ss(zero);
    auto h1_m1_s1 = to_hh_mm_ss(hour + std::chrono::minutes{1} + std::chrono::seconds{1});
    auto h1_11_13p666 = to_hh_mm_ss(std::chrono::seconds{4273} + std::chrono::milliseconds{666});
    assert(std::get<0>(zero_to_hh_mm_ss)==0);
    assert(std::get<1>(zero_to_hh_mm_ss)==0);
    assert(std::get<2>(zero_to_hh_mm_ss)==0.);
    assert(std::get<0>(h1_m1_s1)==1);
    assert(std::get<1>(h1_m1_s1)==1);
    assert(std::get<2>(h1_m1_s1)==1.);
    assert(std::get<0>(h1_11_13p666)==1);
    assert(std::get<1>(h1_11_13p666)==11);
    assert(is_near(std::get<2>(h1_11_13p666),13.666,1.e-12));
}

int main(int argc,char * argv[]) {
    std::cout << "============================== t-utils unit-tests set ==============================\n\n";
    AutoProfilerTester auto_profiler_tester{};
    auto_profiler_tester.test_measurement();
    auto_profiler_tester.test_time_conversions();
    std::cout << "\n\n====================================================================================";
    return 0;
}