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
#include "kadath_bin_bh.hpp"
#include "bns_xcts_solver.hpp"
#include "bns_xcts_setup.hpp"
#include "bns_xcts_regrid.hpp"
#include "Solvers/sequences/parameter_sequence.hpp"
#include<filesystem>

/**
 * \addtogroup BNS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

/**
 * @brief bns_xcts_solution_driver
 *
 * Control computation of BNS from setup config/dat file combination
 * filename is pulled from bconfig which should pair with <filename>.dat
 *
 * Also, KADATH at the time of writing, was strict on not allowing trivial
 * construction of base types (Base_tensor, Scalar, Tensor, Space, etc) which
 * means a driver is required to populate the related Solver class.
 *
 * @tparam config_t Configurator type
 * @param bconfig Configurator object
 * @param outputdir output directory for computed ID
 * @return Success/failure
 */
template<class config_t>
int bns_xcts_solution_driver (config_t& bconfig, std::string outputdir);

/**
 * @brief bns_xcts_driver
 * 
 * Control computation of a single BNS including resolution increase
 * 
 * @tparam Res_t resolution sequence type
 * @tparam config_t Configurator type
 * @param bconfig Configurator object
 * @param resolution Resolution sequence
 * @param outputdir output directory for computed ID
 * @return Success/failure
*/
template<class Res_t, class config_t>
int bns_xcts_driver (config_t & bconfig, 
                          Res_t const & resolution,
                          std::string outputdir);

/**
 * @brief Copy details from minimal or full Config file to construct the full
 * sequence Config file
 * 
 * @tparam config_t Config type
 * @param seqconfig User provided sequence Config
 * @param outputdir location to save files to
 * @return Complete Config file
*/
template<class config_t>
config_t bns_xcts_sequence_setup (config_t & seqconfig, std::string outputdir);

/**
 * @brief Sequence controller for computing binary sequences from scratch
 * 
 * @tparam Seq_t Parameter_sequence type for a binary sequence
 * @tparam Res_t Parameter_sequence type for resolution sequences
 * @param seqconfig User provided sequence Config
 * @param seq 
 * @param outputdir location to save files to
*/
template<class Seq_t, class Res_t, class config_t>
int bns_xcts_sequence (config_t & seqconfig, 
                          Seq_t const & seq,
                          Res_t const & resolution,
                          std::string outputdir);
}}
/** @}*/
#include "bns_xcts_driver.cpp"