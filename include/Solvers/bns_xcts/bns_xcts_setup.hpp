/*
 * Copyright 2023
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
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
#include "mpi.h"
#include "bns_xcts_solver.hpp"
#include <array>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;
/**
 * \addtogroup BNS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

/**
 * bns_xcts_setup_bin
 *
 * Ensure proper binary settings
 *
 * @param[input] bconfig: BNS Configurator file
 */
template<class config_t>
inline void bns_xcts_setup_bin_config(config_t& bconfig);

/**
 * bns_xcts_setup_space
 *
 * Create the numerical space and initial guess for the BNS by
 * obtaining the isolated BH solutions, updating the BNS
 * config file, and importing the BH solutions into the BNS space.
 *
 * @param[input] bconfig: BNS Configurator file
 */
template<class config_t>
void bns_xcts_setup_space (config_t& bconfig);

/**
 * bns_xcts_superimposed_import
 *
 * Flow control to initialize the binary space based on
 * superimposing two isolated compact object solutions.
 * These are simply read from file here before being passed
 * to bns_xcts_setup_boosted_3d for computing the initial
 * guess
 *
 * @param[input] bconfig: binary configurator
 * @param[input] NSfilenames: array of filenames for isolated solutions
 */
template<class config_t>
void bns_xcts_superimposed_import(config_t& bconfig,std::array<std::string, 2> NSfilenames);

/**
 * @brief Generate superimposed guess from 3D isolated solutions
 * 
 * @param NS1config First NS solution
 * @param NS2config Second NS solution
 * @param bconfig Binary Config
 */
template<typename eos_t>
inline void bns_setup_boosted_3d(
  kadath_config_boost<BCO_NS_INFO>& NS1config, kadath_config_boost<BCO_NS_INFO>& NS2config,
  kadath_config_boost<BIN_INFO>& bconfig);
/** @}*/
}}
#include "bns_xcts_setup.cpp"