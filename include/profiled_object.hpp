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
#ifndef __PROFILED_OBJECT_HPP_
#define __PROFILED_OBJECT_HPP_

#include <chrono>
#include <map>
#include <tuple>
#include <typeinfo>
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <deque>

#include <type_traits>
#include <functional>
#include <cmath>
#include "config.h"


namespace std{
    //! Trivial overload for compatibility.
    inline std::string const & to_string(std::string const &s) {return s;}
    //! Trivial overload for compatibility.
    inline std::string to_string(char const * s) {return std::string{s};}
}


namespace Kadath {

    enum : bool {WITHOUT_UNITS = false,WITH_UNITS = true};

    template<typename Duration> struct Time_unit_traits {};
    template<typename T> struct Time_unit_traits<std::chrono::duration<T,std::nano>>
    {
        static constexpr char const * const value = "ns";
        static constexpr char const * const complete_value = "nanosecond(s)";
    };
    template<typename T> struct Time_unit_traits<std::chrono::duration<T,std::micro>>
    {
        static constexpr char const * const value = "µs";
        static constexpr char const * const complete_value = "microsecond(s)";
    };
    template<typename T> struct Time_unit_traits<std::chrono::duration<T,std::milli>>
    {
        static constexpr char const * const value = "ms";
        static constexpr char const * const complete_value = "millisecond(s)";
    };
    template<typename T> struct Time_unit_traits<std::chrono::duration<T>>
    {
        static constexpr char const * const value = "s";
        static constexpr char const * const complete_value = "second(s)";
    };
    template<typename T> struct Time_unit_traits<std::chrono::duration<T,std::ratio<60>>>
    {
        static constexpr char const * const value = "min";
        static constexpr char const * const complete_value = "minute(s)";
    };
    template<typename T> struct Time_unit_traits<std::chrono::duration<T,std::ratio<3600>>>
    {
        static constexpr char const * const value = "h";
        static constexpr char const * const complete_value = "hour(s)";
    };
    template<typename T> struct Time_unit_traits<std::chrono::duration<T,std::ratio<3600 * 24>>>
    {
        static constexpr char const * const value = "d";
        static constexpr char const * const complete_value = "day(s)";
    };
    template<typename T> struct Time_unit_traits<std::chrono::duration<T,std::ratio<3600 * 24 * 365>>>
    {
        static constexpr char const * const value = "y";
        static constexpr char const * const complete_value = "year(s)";
    };

    //! Convert the cmake option value into a static boolean value.
    constexpr bool auto_profiling_enabled =
#ifdef ENABLE_INTERNAL_PROFILER
        true;
#else
        false;
#endif
    //! Dumny type used to pass as undefined type template parameter.
    struct UndefinedType {};

    template<typename OutputDuration> class Profiled_object_base
    {
    public:
        using Hash_key = std::size_t;
        using Duration = std::chrono::high_resolution_clock::duration;

        //! Agglomerate for Statistics.
        struct Statistics {
            std::string user_key;
            unsigned long n_samples;
            double total_duration;
            double average_duration;
            double std_deviation;
        };
        using Stat_map = std::map<Hash_key,Statistics>;

    protected:
        //! Map containing Statistics concerning each sequence of measures.
        static Stat_map statistic_map;

    public:
        //! Read-only access to Statistics.
        static std::map<Hash_key,Statistics> const & get_statistic_map() {return statistic_map;}
        //! Clear all computed statistics.
        static void reset_statistics() {statistic_map.clear();}
        static constexpr char const * get_unit() {return Time_unit_traits<OutputDuration>::value;}
        static constexpr char const * get_complete_unit() {return Time_unit_traits<OutputDuration>::complete_value;}

        /**
         * Sends the data into the passed output stream.
         * @param os output stream to send the data to (can be \c cout, \c cerr, an \c ofstream, etc.).
         * @return a reference toward the current object (\c *this).
         */
        static void display(std::ostream & os,char const * separator = "|",bool with_units = WITH_UNITS);
        /**
         * Helper function for the conversion to different time units.
         * @tparam D output time Duration type (must be an instantiation of the \c std::chrono::Duration template type).
         * @tparam T output numeric representation type (must either be an unsigned integer type or a floating point).
         * @param d time Duration to convert.
         * @return converted time Duration in the desired numeric type.
         */
        template<typename D,typename T=double> static inline T to(Duration const &d)
        {
            static_assert(std::is_floating_point<T>::value || std::is_unsigned<T>::value);
            return std::chrono::duration_cast<std::chrono::duration<T,typename D::period>>(d).count();
        }
        /**
         * Covenience function for the conversion to hours.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time Duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_hours(Duration const &d)
        {return to<std::chrono::hours,T>(d);}
        /**
         * Covenience function for the conversion to minutes.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time Duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_minutes(Duration const &d)
        {return to<std::chrono::minutes,T>(d);}
        /**
         * Covenience function for the conversion to seconds.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time Duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_seconds(Duration const &d)
        {return to<std::chrono::seconds,T>(d);}
        /**
         * Covenience function for the conversion to milliseconds.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time Duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_milliseconds(Duration const &d)
        {return to<std::chrono::milliseconds,T>(d);}
        /**
         * Covenience function for the conversion to microseconds.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time Duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_microseconds(Duration const &d)
        {return to<std::chrono::microseconds,T>(d);}
        /**
         * Covenience function for the conversion to nanoseconds.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time Duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_nanoseconds(Duration const &d)
        {return to<std::chrono::nanoseconds,T>(d);}
        //! Conversion to the hh:mm:ss format.
        static inline std::tuple<unsigned,unsigned,double> to_hh_mm_ss(Duration const & d)
        {
            unsigned const hh {to<std::chrono::hours,unsigned>(d)};
            unsigned const hh_to_min {hh*60};
            unsigned const mm {to<std::chrono::minutes,unsigned>(d) - hh_to_min};
            double const ss {to<std::chrono::seconds>(d) - (hh_to_min + mm)*60.};
            return {hh,mm,ss};
        }
    };

    /**
     * Class used to evaluate the total Duration of multiple calls of some parts of the methods within its derived class,
     * store it in a map with a string keys and perform time units conversions. The profiling operation are deactivated
     * and turned into simple chronometer, without time recording and referencing of the timed pieces of code, when the
     * ENABLED boolean template parameter is false. This profiler class is NOT THREAD SAFE (though it should work with
     * the non-threaded MPI parallel versions of the library).
     * @tparam Derived When defining a class as derived from Profiled_object, pass the type of that class.
     * @tparam ENABLED boolean value used to activate or deactivate the profiling operations.
     */
    template<   class Derived = UndefinedType,
                typename OutputDuration = std::chrono::duration<double>,
                bool ENABLED = auto_profiling_enabled   >
    class Profiled_object;

    //! Specialization for the disabled profiling case.
    template<class Derived,typename OutputDuration> class Profiled_object<Derived,OutputDuration,false> :
            public Profiled_object_base<OutputDuration>
    {
    public:
        using Hash_key = typename Profiled_object_base<OutputDuration>::Hash_key;
        using Clock = std::chrono::high_resolution_clock;
        using Duration = Clock::duration ;
        using Time_point = Clock::time_point;

        static void display(std::ostream &,char const * separator = "|",bool with_unit = WITH_UNITS)  {}
    protected:
        mutable Time_point start;

    public:
        template<typename... T> inline Hash_key start_chrono(T... ) const { start = Clock::now();};
        inline Duration stop_chrono(Hash_key) const {return Clock::now() - start;}
        template<typename... T> inline Duration stop_chrono(T... ) const {return Duration{};};
        void reset() const {}
        void finalize_profiling() const {}
    };

    //! Specialization for the enabled profiling case.
    template<class Derived,typename OutputDuration> class Profiled_object<Derived,OutputDuration,true> :
            public Profiled_object_base<OutputDuration>
    {
    public:
        using Base = Profiled_object_base<OutputDuration>;
        using Hash_key = typename Base::Hash_key ;
        using Statistics = typename Base::Statistics;
        //! Highest resolution STL Clock type.
        using Clock = std::chrono::high_resolution_clock;
        //! The default Duration Clock type.
        using Duration = Clock::duration ;
        //! The default time point type.
        using Time_point = Clock::time_point;
        //! Type storing the computational Duration measures
        using Duration_deque = std::deque<Duration>;
        //! Iterator type for the profiling map.
        using pm_iterator = typename std::map<Hash_key,Duration_deque>::iterator;
        //! Const iterator type for the profiling map.
        using pm_const_iterator = typename std::map<Hash_key,Duration_deque>::const_iterator;
        //! Iterator type for the user key map.
        using uk_iterator = typename std::map<Hash_key,std::string>::iterator;
        //! Const iterator type for the user key map.
        using uk_const_iterator = typename std::map<Hash_key,std::string>::const_iterator;

    protected:
        //! The hash function.
        std::hash<std::string> hash;
        //! Map containing the measured durations, associated to its hash keys.
        mutable std::map<Hash_key,Duration_deque> profiling_map;
        //! Correspondence beetween hash keys and user-given string keys.
        mutable std::map<Hash_key,std::string> user_keys;
        //! Map containing the measures being processed (this design allows overlapping block code measures).
        mutable std::map<Hash_key,Time_point> current_measures;
        //! Name given to the instance of the profiler.
        std::string name;

        //! Variadic overload of the STL to_string function.
        template<typename... T> static std::string to_string(T... args)
        {
            return (std::string{} + ... + std::to_string(args));
        }

    public:
        //! Constructor.
        Profiled_object(char const * _name = typeid(Derived).name()) : hash{}, profiling_map{},
                                                                       user_keys{}, current_measures{}, name{_name} {}
        //! Destructor.
        virtual ~Profiled_object() { this->finalize_profiling();}

        /**
         * Method to call at the begining of a code block to time. The argument are cast into strings and then
         * concatenated to form the key. On the first call, the key is used to insert a new entry in the measures map.
         * @tparam T types of the arguments
         * @param key_prefix
         * @param key_parts
         * @return the hash key calculated from the user-given string.
         */
        template<typename... T> inline Hash_key start_chrono(std::string const & key_prefix, T... key_parts) const
        {
            return _start_chrono(key_prefix + to_string(key_parts...));
        }
        /**
         * Method to call to mark the end of the piece of code to time. The measure time is added to the total time
         * in the measures map.
         * @param key the hash key returned by the call to begin.
         * @return the time measure for this evaluation of the timed block code.
         */
        inline Duration stop_chrono(Hash_key key) const {return _stop_chrono(key);}
        /**
         * Version of stop_chrono using the same arguments as the corresponding call to \c start_chrono. To avoid
         * errors, the best way is to call the hash-key-using overload of this method.
         * @tparam T types of the arguments
         * @param key_prefix
         * @param key_parts
         * @return the time measure for this evaluation of the timed block code.
         */
        template<typename... T> inline Duration stop_chrono(std::string const & key_prefix, T... key_parts) const
        {
            std::string user_key{key_prefix + to_string(key_parts...)};
            return _stop_chrono(hash(name + "." + user_key));
        }

        std::hash<std::string> const & get_hash() const {return hash;}
        //! Accessor for the \c name data member.
        const std::string & get_name() const {return name;}
        //! Mutator for the \c name data member.
        Profiled_object & set_name(const std::string &_name) { name = name; return *this;}
        //! Read-only access to the profiling results map.
        std::map<Hash_key,Duration_deque> const & get_profiling_map() const {return profiling_map;}
        //! Read-only access to the user keys dictionnary.
        std::map<Hash_key,std::string> const & get_user_keys() const {return user_keys;}


    public:
        /**
         * Empty all maps to reset the object.
         * @return Reference toward the current object.
         */
        void reset() const;

         /**
          * Computes statistical data and erase measure-related maps. Another set of measure can be made, and if so,
          * another call to \c finalize_profiling() will update the Statistics according to the new data.
          * @return a reference toward the current object.
          */
         void finalize_profiling() const;

    private:
        //! Implementation of \c start_chrono (see the interface without the underscore prefix for informations).
        Hash_key _start_chrono(std::string const & user_key) const;
        //! Implementation of \c stop_chrono (see the interface without the underscore prefix for informations).
        Duration _stop_chrono(Hash_key key) const;

    };


    template<class Derived,typename OD>
    typename Profiled_object<Derived,OD,true>::Hash_key
    Profiled_object<Derived,OD,true>::_start_chrono(std::string const &user_key) const
    {
        Hash_key const key {hash(name + "." + user_key)};
        auto user_key_it = user_keys.find(key);
        pm_iterator pm_location{};
        if(user_key_it == user_keys.end())
        {
            user_keys[key] = user_key;
            auto insertion_result = profiling_map.emplace(key, Duration_deque{});
            assert(insertion_result.second);
            pm_location = insertion_result.first;
        } else
        {
            pm_location = profiling_map.find(key);
            assert(user_key_it->second == user_key);
        }
        Duration_deque & measure = pm_location->second;

        assert(current_measures.find(key) == current_measures.end());
        if(current_measures.find(key) != current_measures.end())
        {
            std::cerr << "WARNING: while profiling " << name
                      << " a measure with the same key already exists. This may be due either by different threads trying to"
                         "measure at the same time, or overlapping block with the same key. In the later case, try changing "
                         "one of the key." << std::endl;
            std::cerr << " ==> Problematic measure : " << user_key << std::endl;
        }
        Time_point & id_current_measures = current_measures[key];
        id_current_measures = Clock::now();
        return key;
    }

    template<class Derived,typename OD>
    typename Profiled_object<Derived,OD,true>::Duration Profiled_object<Derived,OD,true>::_stop_chrono(Hash_key key) const
    {
        Time_point const current_measure_stop_time {Clock::now()};
        typename std::map<Hash_key,Time_point>::iterator const current_measure{current_measures.find(key)};
        if(current_measure == current_measures.end())
        {
            std::cerr << "Error while profiling " << name << " measure key cannot be "
                                                             "found (maybe start_chrono() has been forgotten or called before stop_chrono()). "
                                                             "The current measure will be discarded." << std::endl;
            std::cerr << " ==> Problematic measure : " << user_keys[key] << std::endl;
            return Duration{};
        }
        else
        {
            Duration const current_measure_duration {current_measure_stop_time - current_measure->second};
            profiling_map[key].push_back(current_measure_duration);
            current_measures.erase(current_measure);
            return current_measure_duration;
        }
    }

    template<class Derived,typename OD> void Profiled_object<Derived,OD,true>::reset() const
    {
        profiling_map.clear();
        user_keys.clear();
        current_measures.clear();
    }

    template<class Derived,typename OD> void Profiled_object<Derived,OD,true>::finalize_profiling() const
    {
        assert(current_measures.empty());
        assert(profiling_map.size() == user_keys.size());
        pm_const_iterator i {profiling_map.begin()};
        uk_const_iterator j {user_keys.begin()};
        for(;i != profiling_map.end() && j != user_keys.end();i++,j++)
        {
            Hash_key const key {i->first};
            assert(key == j->first);
            typename std::map<Hash_key,Statistics>::iterator stat_it = Base::statistic_map.find(key);
            if(stat_it == Base::statistic_map.end())
            {
                std::string input_name {name};
                input_name += ".";
                input_name += j->second;
                std::tie(stat_it,std::ignore) =
                    Base::statistic_map.emplace(key, Statistics{input_name, 0ul, 0., 0., 0.});
            }
            unsigned long const n1{stat_it->second.n_samples}, n2{static_cast<unsigned long>(i->second.size())};
            double const total1 {stat_it->second.total_duration};
            double const mu1 {stat_it->second.average_duration};
            double const total2 {std::accumulate(i->second.begin(),i->second.end(),0.,
                                        [](double l,Duration const & r){return l + Base::template to<OD>(r);})};
            double const mu2 {total2 / n2};
            double const sigma1_2 {stat_it->second.std_deviation};
            double const n2_sigma2_2 {std::accumulate(i->second.begin(),i->second.end(),0.,
                                            [mu2](double l,Duration const &r) {
                                                double const d_r_mu2{Base::template to<OD>(r) - mu2};
                                                return l + d_r_mu2 * d_r_mu2;})};
            unsigned long const ntot {n1 + n2};
            stat_it->second.n_samples = ntot;
            stat_it->second.total_duration = total1 + total2;
            stat_it->second.average_duration = (n1*mu1 + n2*mu2) / ntot;
            stat_it->second.std_deviation = (n1*sigma1_2 + n2_sigma2_2 + n1*n2*(mu1-mu2)*(mu1-mu2)/ntot) / ntot;
        }
        this->reset();
    }

    template<typename OD>
    void Profiled_object_base<OD>::display(std::ostream &os, char const *separator, bool with_units)
    {
        // entries are sorted from the greater total Duration to the smallest.
        std::vector<std::pair<Hash_key,Statistics>> sorted_values(statistic_map.begin(), statistic_map.end());
        std::sort(sorted_values.begin(),sorted_values.end(),[](std::pair<Hash_key,Statistics> const &l,
                std::pair<Hash_key,Statistics> const &r) {return l.second.total_duration > r.second.total_duration;});
        auto  insert_unit_and_separator = [separator,with_units,&os](bool sep = true)
                {
                    if( with_units)
                    {
                        os << ' ' << std::left << std::setw(3) << get_unit();
                    } else {
                        os << "    ";
                    }
                    if(sep) os << ' ' << separator << ' ';
                };
        std::deque<std::pair<Hash_key,std::string>> too_long_id;
        for (std::pair<Hash_key, Statistics> const &data : sorted_values)
        {
            std::string const id {data.second.user_key};
            if(id.size() < 60) {
                os << std::right << std::setw(60) << data.second.user_key;
            } else {
                Hash_key const id_hk {std::hash<std::string>{}(id)};
                too_long_id.push_back({id_hk,id});
                os << std::right << std::setw(60) << id_hk;
            }
            os << " " << separator << ' ';
            os << std::fixed << std::setprecision(9) << std::setw(10) << data.second.total_duration;
            insert_unit_and_separator();
            os << std::fixed << std::setprecision(9) << std::setw(10) << data.second.average_duration;
            insert_unit_and_separator();
            os << std::fixed << std::setprecision(9) << std::setw(10) << sqrt(data.second.std_deviation);
            insert_unit_and_separator();
            os << std::fixed << std::setprecision(9) << std::setw(10) << data.second.n_samples << ' ';
            if(with_units) os << "sample(s) ";
            os << "\n";
        }
        for(auto const & tlid : too_long_id)
        {
            os << "#" << std::setw(20) << tlid.first << " : " << tlid.second  << std::endl;
        }
    }

    template<class C> void profiling_report(C const & c,std::ostream & os)
    {
#ifdef ENABLE_INTERNAL_PROFILER
        os << "=======================================================================" << std::endl;
        os << "Profiling report : " << std::endl;
        os << "=======================================================================" << std::endl << std::endl;
        os << "                            id                              " << " | ";
        os << "  total time    |  average  time  |    std. dev.    |    N " << std::endl;
        os << "============================================================" << "=|=";
        os << "================|=================|=================|======================" << std::endl;
#endif
        c.display(std::cout);
    }


}

#endif //__PROFILED_OBJECT_HPP_
