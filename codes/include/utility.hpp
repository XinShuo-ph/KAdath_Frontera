//
// Created by sauliac on 29/05/2020.
//
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

#ifndef __KADATH_CODES_UTILITY_HPP_
#define __KADATH_CODES_UTILITY_HPP_

/****************************************************************************************
 *  Some usefull preprocessor macros for class with a lot of data members.
 ****************************************************************************************/

#define decl_data_member_with_accessor_and_mutator_by_val(type,identifier) \
protected: \
    type identifier; \
public: \
    type get_##identifier () const {return identifier;} \
    void set_##identifier(type new_value) { \
        identifier = new_value; \
    }

#define decl_data_member_with_accessor_by_val(type,identifier) \
protected: \
    type identifier; \
public: \
    type get_##identifier () const {return identifier;}


#define decl_data_member_with_accessor_and_mutator_by_ref(type,identifier) \
protected: \
    type identifier; \
public: \
    type const & get_##identifier () const {return identifier;} \
    void set_##identifier(type const & new_value) { \
        identifier = new_value; \
    }

#define decl_data_member_with_accessor_by_ref(type,identifier) \
protected: \
    type identifier; \
public: \
    type const & get_##identifier () const {return identifier;}

#define decl_ptr_member_with_rw_accessor(type,identifier) \
protected: \
    std::unique_ptr<type> identifier; \
public: \
    std::unique_ptr<type> const & get_##identifier () const {return identifier;} \
    std::unique_ptr<type> & get_##identifier () {return identifier;}

#define decl_ptr_member_with_r_accessor(type,identifier) \
protected: \
    std::unique_ptr<type> identifier; \
public: \
    std::unique_ptr<type> const & get_##identifier () const {return identifier;}



/****************************************************************************************
 * A simple argument parser.
 ****************************************************************************************/

template<typename T> struct is_admissible_option_type {
    static constexpr bool value = false;
};
template<> struct is_admissible_option_type<bool> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<int> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<long> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<long long> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<double> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<std::string> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<char*> {static constexpr bool value = true;};

template<typename T> inline T from_string(std::string const & s) {
    std::stringstream sss{s};
    T t_from_s{};
    sss >> t_from_s;
    return t_from_s;
}
template<> inline std::string from_string<std::string>(std::string const & s) {return s;}
template<> inline char const * from_string<char const *>(std::string const & s) {return s.c_str();}


class Arguments_parser{
private:
    std::vector<std::string> command_line;

public:
    Arguments_parser (int &argc, char **argv){
        for (int i=1; i < argc; ++i)
            this->command_line.push_back(std::string(argv[i]));
    }

    template<typename OT>
    std::pair<typename std::enable_if<
                is_admissible_option_type<OT>::value,
                OT
            >::type,bool>
    get_option_value(const std::string &option) const {
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->command_line.begin(), this->command_line.end(), option);
        if (itr != this->command_line.end() && ++itr != this->command_line.end()){
            return {from_string<OT>(*itr),true};
        }
        else return {OT{},false};
    }
};



#endif //__KADATH_CODES_UTILITY_HPP_
