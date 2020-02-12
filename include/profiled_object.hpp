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
#include <cassert>
#include <type_traits>
#include <functional>
#include "config.h"


namespace std{
    //! Trivial overload for compatibility.
    inline std::string const & to_string(std::string const &s) {return s;}
    //! Trivial overload for compatibility.
    inline std::string to_string(char const * s) {return std::string{s};}
}


namespace Kadath {

    enum : bool {WITHOUT_UNITS = false,WITH_UNITS = true};

    template<typename Duration> struct TimeUnitTraits {};
    template<typename T> struct TimeUnitTraits<std::chrono::duration<T,std::nano>>
    {
        static constexpr char const * const value = "ns";
        static constexpr char const * const complete_value = "nanosecond(s)";
    };
    template<typename T> struct TimeUnitTraits<std::chrono::duration<T,std::micro>>
    {
        static constexpr char const * const value = "µs";
        static constexpr char const * const complete_value = "microsecond(s)";
    };
    template<typename T> struct TimeUnitTraits<std::chrono::duration<T,std::milli>>
    {
        static constexpr char const * const value = "ms";
        static constexpr char const * const complete_value = "millisecond(s)";
    };
    template<typename T> struct TimeUnitTraits<std::chrono::duration<T>>
    {
        static constexpr char const * const value = "s";
        static constexpr char const * const complete_value = "second(s)";
    };
    template<typename T> struct TimeUnitTraits<std::chrono::duration<T,std::ratio<60>>>
    {
        static constexpr char const * const value = "min";
        static constexpr char const * const complete_value = "minute(s)";
    };
    template<typename T> struct TimeUnitTraits<std::chrono::duration<T,std::ratio<3600>>>
    {
        static constexpr char const * const value = "h";
        static constexpr char const * const complete_value = "hour(s)";
    };
    template<typename T> struct TimeUnitTraits<std::chrono::duration<T,std::ratio<3600*24>>>
    {
        static constexpr char const * const value = "d";
        static constexpr char const * const complete_value = "day(s)";
    };
    template<typename T> struct TimeUnitTraits<std::chrono::duration<T,std::ratio<3600*24*365>>>
    {
        static constexpr char const * const value = "y";
        static constexpr char const * const complete_value = "year(s)";
    };

    //! Convert the cmake option value into a static boolean value.
    constexpr bool auto_profiling_enabled {ENABLE_PROFILING == 1};
    //! Dumny type used to pass as undefined type template parameter.
    struct UndefinedType {};

    template<typename OutputDuration> class ProfiledObjectBase
    {
    public:
        using hash_key = std::size_t;
        using duration = std::chrono::high_resolution_clock::duration;

        //! Agglomerate for statistics.
        struct statistics {
            std::string user_key;
            unsigned long n_samples;
            double total_duration;
            double average_duration;
            double std_deviation;
        };

    protected:
        //! Map containing statistics concerning each sequence of measures.
        inline static std::map<hash_key,statistics> statistic_map{};

    public:
        //! Read-only access to statistics.
        static std::map<hash_key,statistics> const & get_statistic_map() {return statistic_map;}

        static constexpr char const * get_unit() {return TimeUnitTraits<OutputDuration>::value;}
        static constexpr char const * get_complete_unit() {return TimeUnitTraits<OutputDuration>::complete_value;}

        /**
         * Sends the data into the passed output stream.
         * @param os output stream to send the data to (can be \c cout, \c cerr, an \c ofstream, etc.).
         * @return a reference toward the current object (\c *this).
         */
        static void display(std::ostream & os,char const * separator = "|",bool with_units = WITH_UNITS);
        /**
         * Helper function for the conversion to different time units.
         * @tparam D output time duration type (must be an instantiation of the \c std::chrono::duration template type).
         * @tparam T output numeric representation type (must either be an unsigned integer type or a floating point).
         * @param d time duration to convert.
         * @return converted time duration in the desired numeric type.
         */
        template<typename D,typename T=double> static inline T to(duration const &d)
        {
            static_assert(std::is_floating_point<T>::value || std::is_unsigned<T>::value);
            return std::chrono::duration_cast<std::chrono::duration<T,typename D::period>>(d).count();
        }
        /**
         * Covenience function for the conversion to hours.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_hours(duration const &d)
        {return to<std::chrono::hours,T>(d);}
        /**
         * Covenience function for the conversion to minutes.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_minutes(duration const &d)
        {return to<std::chrono::minutes,T>(d);}
        /**
         * Covenience function for the conversion to seconds.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_seconds(duration const &d)
        {return to<std::chrono::seconds,T>(d);}
        /**
         * Covenience function for the conversion to milliseconds.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_milliseconds(duration const &d)
        {return to<std::chrono::milliseconds,T>(d);}
        /**
         * Covenience function for the conversion to microseconds.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_microseconds(duration const &d)
        {return to<std::chrono::microseconds,T>(d);}
        /**
         * Covenience function for the conversion to nanoseconds.
         * @tparam T numeric type to express the conversion result to (\c double by default).
         * @param d the time duration to convert.
         * @return the result of the conversion in the desired numeric type.
         */
        template<typename T=double> static inline T to_nanoseconds(duration const &d)
        {return to<std::chrono::nanoseconds,T>(d);}
        //! Conversion to the hh:mm:ss format.
        static inline std::tuple<unsigned,unsigned,double> to_hh_mm_ss(duration const & d)
        {
            unsigned const hh {to<std::chrono::hours,unsigned>(d)};
            unsigned const hh_to_min {hh*60};
            unsigned const mm {to<std::chrono::minutes,unsigned>(d) - hh_to_min};
            double const ss {to<std::chrono::seconds>(d) - (hh_to_min + mm)*60.};
            return {hh,mm,ss};
        }
    };

    /**
     * Class used to evaluate the total duration of multiple calls of some parts of the methods within its derived class,
     * store it in a map with a string keys and perform time units conversions. The profiling operation are deactivated
     * and turned into zero-cost do-nothing methods when the ENABLED boolean template parameter is false. Keep in mind
     * this profiler class is NOT THREAD SAFE (though it should work with the non-threaded MPI parallel versions of
     * the library).
     * @tparam Derived When defining a class as derived from ProfiledObject, pass the type of that class.
     * @tparam ENABLED boolean value used to activate or deacivate the profiling operations.
     */
    template<   class Derived = UndefinedType,
                typename OutputDuration = std::chrono::duration<double>,
                bool ENABLED = auto_profiling_enabled   >
    class ProfiledObject;

    //! Specialization for the disabled profiling case.
    template<class Derived,typename OutputDuration> class ProfiledObject<Derived,OutputDuration,false> :
            public ProfiledObjectBase<OutputDuration>
    {
    public:
        using hash_key = std::size_t;
        using clock = std::chrono::high_resolution_clock;
        using duration = clock::duration ;
        using time_point = clock::time_point;

        static void display(std::ostream &,char const * separator = "|",bool with_unit = WITH_UNITS)  {}
    protected:
        mutable time_point start;

    public:
        template<typename... T> inline hash_key start_chrono(T... ) const {start = clock::now();};
        inline duration stop_chrono(hash_key) const {return clock::now()-start;}
        template<typename... T> inline duration stop_chrono(T... ) const {return duration{};};
        void reset() const {}
        void finalize() const {}
    };

    //! Specialization for the enabled profiling case.
    template<class Derived,typename OutputDuration> class ProfiledObject<Derived,OutputDuration,true> :
            public ProfiledObjectBase<OutputDuration>
    {
    public:
        using Base = ProfiledObjectBase<OutputDuration>;
        using hash_key = typename Base::hash_key ;
        using statistics = typename Base::statistics;
        //! Highest resolution STL clock type.
        using clock = std::chrono::high_resolution_clock;
        //! The default duration clock type.
        using duration = clock::duration ;
        //! The default time point type.
        using time_point = clock::time_point;
        //! Type storing the computational duration measures
        using duration_deque = std::deque<duration>;
        //! Iterator type for the profiling map.
        using pm_iterator = typename std::map<hash_key,duration_deque>::iterator;
        //! Const iterator type for the profiling map.
        using pm_const_iterator = typename std::map<hash_key,duration_deque>::const_iterator;
        //! Iterator type for the user key map.
        using uk_iterator = typename std::map<hash_key,std::string>::iterator;
        //! Const iterator type for the user key map.
        using uk_const_iterator = typename std::map<hash_key,std::string>::const_iterator;

    protected:
        //! The hash function.
        std::hash<std::string> hash;
        //! Map containing the measured durations, associated to its hash keys.
        mutable std::map<hash_key,duration_deque> profiling_map;
        //! Correspondence beetween hash keys and user-given string keys.
        mutable std::map<hash_key,std::string> user_keys;
        //! Map containing the measures being processed (this design allows overlapping block code measures).
        mutable std::map<hash_key,time_point> current_measures;
        //! Name given to the instance of the profiler.
        std::string name;

        //! Variadic overload of the STL to_string function.
        template<typename... T> static std::string to_string(T... args)
        {
            return (std::string{} + ... + std::to_string(args));
        }

    public:
        //! Constructor.
        ProfiledObject(char const * _name = typeid(Derived).name()) : hash{}, profiling_map{},
            user_keys{}, current_measures{}, name{_name} {}
        //! Destructor.
        virtual ~ProfiledObject() { this->finalize();}

        /**
         * Method to call at the begining of a code block to time. The argument are cast into strings and then
         * concatenated to form the key. On the first call, the key is used to insert a new entry in the measures map.
         * @tparam T types of the arguments
         * @param key_prefix
         * @param key_parts
         * @return the hash key calculated from the user-given string.
         */
        template<typename... T> inline hash_key start_chrono(std::string const & key_prefix, T... key_parts) const
        {
            return _start_chrono(key_prefix + to_string(key_parts...));
        }
        /**
         * Method to call to mark the end of the piece of code to time. The measure time is added to the total time
         * in the measures map.
         * @param key the hash key returned by the call to begin.
         * @return the time measure for this evaluation of the timed block code.
         */
        inline duration stop_chrono(hash_key key) const {return _stop_chrono(key);}
        /**
         * Version of stop_chrono using the same arguments as the corresponding call to \c start_chrono. To avoid
         * errors, the best way is to call the hash-key-using overload of this method.
         * @tparam T types of the arguments
         * @param key_prefix
         * @param key_parts
         * @return the time measure for this evaluation of the timed block code.
         */
        template<typename... T> inline duration stop_chrono(std::string const & key_prefix, T... key_parts) const
        {
            std::string user_key{key_prefix + to_string(key_parts...)};
            return _stop_chrono(hash(name + "." + user_key));
        }

        //! Accessor for the \c name data member.
        const std::string & get_name() const {return name;}
        //! Mutator for the \c name data member.
        ProfiledObject & set_name(const std::string &_name) {name = name; return *this;}
        //! Read-only access to the profiling results map.
        std::map<hash_key,duration_deque> const & get_profiling_map() const {return profiling_map;}
        //! Read-only access to the user keys dictionnary.
        std::map<hash_key,std::string> const & get_user_keys() const {return user_keys;}





    public:
        /**
         * Empty all maps to reset the object.
         * @return Reference toward the current object.
         */
        void reset() const;

         /**
          * Computes statistical data and erase measure-related maps. Another set of measure can be made, and if so,
          * another call to \c finalize() will update the statistics according to the new data.
          * @return a reference toward the current object.
          */
         void finalize() const;

    private:
        //! Implementation of \c start_chrono (see the interface without the underscore prefix for informations).
        hash_key _start_chrono(std::string const & user_key) const;
        //! Implementation of \c stop_chrono (see the interface without the underscore prefix for informations).
        duration _stop_chrono(hash_key key) const;

    };


    template<class Derived,typename OD>
    typename ProfiledObject<Derived,OD,true>::hash_key
    ProfiledObject<Derived,OD,true>::_start_chrono(std::string const &user_key) const
    {
        hash_key const key {hash(name + "." + user_key)};
        auto user_key_it = user_keys.find(key);
        pm_iterator pm_location{};
        if(user_key_it == user_keys.end())
        {
            user_keys[key] = user_key;
            auto insertion_result = profiling_map.emplace(key, duration_deque{});
            assert(insertion_result.second);
            pm_location = insertion_result.first;
        } else
        {
            pm_location = profiling_map.find(key);
            assert(user_key_it->second == user_key);
        }
        duration_deque & measure = pm_location->second;

        assert(current_measures.find(key) == current_measures.end());
        if(current_measures.find(key) != current_measures.end())
        {
            std::cerr << "WARNING: while profiling " << name
                      << " a measure with the same key already exists. This may be due either by different threads trying to"
                         "measure at the same time, or overlapping block with the same key. In the later case, try changing "
                         "one of the key." << std::endl;
            std::cerr << " ==> Problematic measure : " << user_key << std::endl;
        }
        time_point & id_current_measures = current_measures[key];
        id_current_measures = clock::now();
        return key;
    }

    template<class Derived,typename OD>
    typename ProfiledObject<Derived,OD,true>::duration ProfiledObject<Derived,OD,true>::_stop_chrono(hash_key key) const
    {
        time_point const current_measure_stop_time {clock::now()};
        typename std::map<hash_key,time_point>::iterator const current_measure{current_measures.find(key)};
        if(current_measure == current_measures.end())
        {
            std::cerr << "Error while profiling " << name << " measure key cannot be "
                                                             "found (maybe start_chrono() has been forgotten or called before stop_chrono()). "
                                                             "The current measure will be discarded." << std::endl;
            std::cerr << " ==> Problematic measure : " << user_keys[key] << std::endl;
            return duration{};
        }
        else
        {
            duration const current_measure_duration {current_measure_stop_time - current_measure->second};
            profiling_map[key].push_back(current_measure_duration);
            current_measures.erase(current_measure);
            return current_measure_duration;
        }
    }

    template<class Derived,typename OD> void ProfiledObject<Derived,OD,true>::reset() const
    {
        profiling_map.clear();
        user_keys.clear();
        current_measures.clear();
    }

    template<class Derived,typename OD> void ProfiledObject<Derived,OD,true>::finalize() const
    {
        assert(current_measures.empty());
        assert(profiling_map.size() == user_keys.size());
        pm_const_iterator i {profiling_map.begin()};
        uk_const_iterator j {user_keys.begin()};
        for(;i != profiling_map.end() && j != user_keys.end();i++,j++)
        {
            hash_key const key {i->first};
            assert(key == j->first);
            typename std::map<hash_key,statistics>::iterator stat_it = Base::statistic_map.find(key);
            if(stat_it == Base::statistic_map.end())
            {
                std::string input_name {name};
                input_name += ".";
                input_name += j->second;
                std::tie(stat_it,std::ignore) =
                    Base::statistic_map.emplace(key,statistics{input_name,0ul,0.,0.,0.});
            }
            unsigned long nsamples {stat_it->second.n_samples + static_cast<unsigned long>(i->second.size())};
            double total {stat_it->second.total_duration};
            double average {total};
            double std_dev {stat_it->second.std_deviation * stat_it->second.n_samples};
            total = std::accumulate(i->second.begin(),i->second.end(),total,[](double l,duration const & r)
                {return l + Base::template to<OD>(r);});
            average = total;
            average /= nsamples;
            std_dev = std::accumulate(i->second.begin(),i->second.end(),std_dev,[average](double l,duration const &r)
                {double const _r{Base::template to<OD>(r)-average}; return l + _r*_r;});
            std_dev /= nsamples;
            stat_it->second.n_samples = nsamples;
            stat_it->second.total_duration = total;
            stat_it->second.average_duration = average;
            stat_it->second.std_deviation = std_dev;
        }
        this->reset();
    }

    template<typename OD>
    void ProfiledObjectBase<OD>::display(std::ostream &os, char const *separator, bool with_units)
    {
        // entries are sorted from the greater total duration to the smallest.
        std::vector<std::pair<hash_key,statistics>> sorted_values(statistic_map.begin(),statistic_map.end());
        std::sort(sorted_values.begin(),sorted_values.end(),[](std::pair<hash_key,statistics> const &l,
                std::pair<hash_key,statistics> const &r) {return l.second.total_duration > r.second.total_duration;});
        auto  insert_unit_and_separator = [separator,with_units,&os](bool sep = true)
                {
                    if( with_units)
                    {
                        os << ' ' << std::setw(3) << get_unit();
                    } else {
                        os << "    ";
                    }
                    if(sep) os << ' ' << separator << ' ';
                };
        std::deque<std::pair<hash_key,std::string>> too_long_id;
        for (std::pair<hash_key, statistics> const &data : sorted_values)
        {
            std::string const id {data.second.user_key};
            if(id.size() < 60) {
                os << std::setw(60) << data.second.user_key;
            } else {
                hash_key const id_hk {std::hash<std::string>{}(id)};
                too_long_id.push_back({id_hk,id});
                os << std::setw(60) << id_hk;
            }
            os << " " << separator << ' ';
            os << std::setw(10) << data.second.total_duration;
            insert_unit_and_separator();
            os << std::setw(10) << data.second.average_duration;
            insert_unit_and_separator();
            os << std::setw(10) << sqrt(data.second.std_deviation);
            insert_unit_and_separator();
            os << std::setw(10) << data.second.n_samples << ' ';
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
#if ENABLE_PROFILING==1
        os << "=======================================================================" << std::endl;
        if(c.get_statistic_map().empty())
        {
            c.finalize();
        }
        os << "Profiling report : " << std::endl;
        os << "=======================================================================" << std::endl << std::endl;
        os << "                            id                              " << " | ";
        os << " total time    |  average time  |   std. dev.    |    N " << std::endl;
        os << "============================================================" << "=|=";
        os << "===============|================|================|======================" << std::endl;
#endif
        c.display(std::cout);
    }

}

#endif //__PROFILED_OBJECT_HPP_
