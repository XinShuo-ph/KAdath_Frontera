/*
 * This file is part of the KADATH library.
 * Author: Samuel Tootle
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/optional.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include "config_enums.hpp"
#include "config_utils_boost.hpp"
#include <stdexcept>
#include <variant>
#include <iostream>
#include <iomanip>
#include <string>
#include <typeinfo>
#include <sstream>
#include <cmath>
#include <utility>
namespace pt = boost::property_tree;
using Tree = pt::ptree;

template<typename ParamC> class kadath_config_boost;
template<typename ParamC> std::ostream& operator<<(std::ostream& out, const kadath_config_boost<ParamC>& config);
/**
 * \addtogroup Configurator
 * @{*/

/**
 * class kadath_config_boost
 * interface to reading/writing parameters from .info file
 *
 * @tparam ParamC type of parameter container (e.g. BH_INFO, BIN_INFO)
 */
template <typename ParamC>
class kadath_config_boost {
  using FArray = std::array<bool, NUM_BCO_FIELDS>;
  using SArray = std::array<bool, NUM_STAGES>;
  using CArray = std::array<bool, NUM_CONTROLS>;
  private:
    std::string filename{};
    std::string outputdir{"./"};
    Tree tree;     
    Tree branch_in;
    FArray fields{};       
    SArray stages{};       
    CArray controls{};     
    std::array<double, NUM_SEQ_SETTINGS> seq_settings{};
    ParamC container;      ///< Parameter container

  public:
    kadath_config_boost();
    
    /**
      * Default constructor used for reading in *.info config files
      *
      * @tparam ParamC type of parameter container (e.g. BH_INFO, BIN_INFO)
      * @param[input] ifile: input filename
      *
     */
    kadath_config_boost(std::string ifile);
    
    /**
      * kadath_config_boost::open_config
      *
      * read-in configuration from filename
      * this will populate the parameter container.
      *
      * @tparam ParamC type of parameter container (e.g. BH_INFO, BIN_INFO)
      * @return status
     */
    int open_config();
    
    /**
      * kadath_config_boost::write_config()
      * write the current configuration to file
      * Default: write to this->outputdir/this->filename
      * Secondary: write to this->outputdir/ofile
      * Last: write to ofile IFF ofile contains a path (i.e. a '/')
      *
      * @tparam ParamC type of parameter container (e.g. BCO_BH_INFO, BIN_INFO)
      * @param [input] ofile: output filename to use.
     */
    void write_config(std::string ofile="null");

    /**
      * kadath_config_boost::config_filename()
      * @return filename: returns constant filename
      */
    const std::string config_filename() const { return filename; }
    const std::string config_filename_abs() const;                  ///< config filename including path
    const std::string space_filename() const;                      ///< space filename including path
    const std::string config_outputdir()const { return outputdir; }

    FArray const& return_fields() const { return fields; }
    FArray& return_fields() { return fields; }

    SArray const& return_stages() const { return stages; }
    SArray& return_stages() { return stages; }

    auto& set_field(const int idx) { return fields[idx]; }
    auto const & field(const int idx) const { return fields[idx]; }

    auto& set_stage(const int idx) { return stages[idx]; }
    auto const & set_stage(const int idx) const { return stages[idx]; }

    auto const & control(const int idx) const { return controls[idx]; }   
    auto& control(const int idx) { return controls[idx]; } 

    auto& seq_setting(const int idx) { return seq_settings[idx]; }
    auto const & seq_setting(const int idx) const { return seq_settings[idx]; }

    void set_outputdir(std::string dir);
    void set_filename(std::string fname);

    /**
      * kadath_config_boost::initialize_binary
      * initialize binary parameter container based on node types
      * @param[input] bco_types: array containing two strings such as {"bh", "bh"}
      */
    void initialize_binary(std::array<std::string, 2> bco_types) { container.init_binary(bco_types); }

    /**
      * kadath_config_boost::set_defaults
      * for a given parameter container, use the conatainers defaults.
      * Note: this is currently used only for isolated objects such as NS or BH.  Not binaries.
      */
    void set_defaults() { container.set_defaults(*this); }

    /**
      * kadath_config_boost::get_map
      * get map of base parameter container
      */
    constexpr auto& get_map() const { 
      return container.get_map();
    }

    /**
      * kadath_config_boost::get_map
      * get map of child parameter container where idx is the child idx (i.e. BCO1)
      * @param[input] idx: child container array index
      */
    constexpr auto& get_map(const int idx) const { 
      return container.get_map(idx);
    }

    /**
      * overloaded () to obtain parameters from base parameter container
      * but throws an error when a NaN is encountered
      * @param[input] idx: index of parameter to read/write
      */
    constexpr auto& operator()(const int idx) { 
      if(!check_for_nan(container.get_map(), container(idx), idx)) {
        return container(idx);
      }
      throw std::invalid_argument("\nInvalid Configuration Indicies\n)");
    }
    /**
      * overloaded () to obtain parameters from child paramter container
      * but throws an error when a NaN is encountered
      * @param[input] idx: index of parameter to read/write
      * @param[input] Pidx: index of child parameter container i.e. [BCO1,BCO2]
      * @return parameter
      */
    constexpr auto& operator()(const int idx, const int Pidx) {
      if(!check_for_nan(container.get_map(Pidx), container(idx, Pidx), idx)) {
        return container(idx, Pidx);
      }
      throw std::invalid_argument("\nInvalid Configuration Indicies\n)");
    }

    /**
      * kadath_config_boost::set()
      * set funct to set parameters from base or child parameter container
      * @tparam idx_t parameter pack type of indexes
      * @param[input] idxs: index/indices of parameter to set
      * @return parameter
      */
    template<typename... idx_t>
    constexpr auto& set(const idx_t... idxs) { 
      return container(idxs...);
    }
   
    /**
      * kadath_config_boost::set_eos()
      * set funct to set parameters of the EOS for base or child parameter container
      * @tparam idx_t parameter pack type of indexes
      * @param[input] idxs: index/indices of eos parameter to set
     */
    template<typename... idx_t>
    constexpr auto& set_eos(const idx_t... idxs) { 
      return container.set_eos_param(idxs...); 
    }

    /**
      * kadath_config_boost::eos()
      * get parameter value of the std::variant eos_parameter for base or child parameter container
      * (see config_bco.hpp for expected types) 
      *
      * here we use template 'specializations' based on type traits to determine
      * how to extract data from a std::variant container.
      * @tparam T type of called parameter
      * @tparam idx_t type of index parameter pack
      * @param[input] idxs: index/indices of eos parameter to read
     */
    template<typename T, typename... idx_t,
    std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
    constexpr const T eos(idx_t... idx) {
      T var{};
      auto set_var = [&](auto&& v) mutable {
        using v_t = std::decay_t<decltype(v)>;  
        if constexpr (std::is_arithmetic_v<v_t>) {
          if(!std::isnan(v))
            var = v;
        } 
      };  
      std::visit(set_var, this->set_eos(idx...));
      return var;
    }

    template<typename T, typename... idx_t,
    std::enable_if_t<std::is_same_v<T, std::string>, bool> = true>
    constexpr const T eos(idx_t... idx) {
      T var{};
      auto set_var = [&](auto&& v) mutable {
        using v_t = std::decay_t<decltype(v)>;  
        if constexpr (std::is_same_v<v_t, std::string>) {   
            var = v;    
        }    
      };  
      std::visit(set_var, this->set_eos(idx...));
      return var;
    }
    
    inline void set_seq_defaults();

    friend std::ostream& operator<< <>(std::ostream&, const kadath_config_boost<ParamC>&) ;
};
/**
 * @}*/

#include "configurator_boost_imp.cpp"
