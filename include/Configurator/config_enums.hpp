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
#include <map>
#include <string>
/**
  * @defgroup Configurator_enums
	* @ingroup Configurator
  * Physically motivated, static enums that allow for flexible reading
  * of boost::ptree while allowing the ability to add/modify/delete params at a
  * later date without having to modify BIN/BCO classes that read in data. This
  * also allows for calling container locations using physically meaningful names
  * @{
  */

/** @brief enum BIN_PARAMS enumerator over binary parameters */
enum BIN_PARAMS {
  DIST=0, DDIST, BIN_RES, GOMEGA, COM, COMY, \
  REXT, QPIG, Q, ADOT, ECC_OMEGA, OUTER_SHELLS, NUM_BPARAMS \
};
/** @brief enum BCO_PARAMS enumerator over BCO parameters */
enum BCO_PARAMS {
  RIN=0, BCO_RES, BCO_QPIG, RMID, FIXED_R, ROUT, NSHELLS, OMEGA, CHI, \
  MIRR, MCH, MB, NC, HC, FIXED_LAPSE, MADM, \
  QLMADM, DIM, USE_TOV1D, FIXED_BCOMEGA, BVELX, BVELY, 
  DECAY, KERR_CHI, KERR_MCH, NINSHELLS, NUM_BCO_PARAMS
};
/** @brief enum EOS_PARAMS enumerator over EOS parameters */
enum EOS_PARAMS {EOSTYPE, EOSFILE, HCUT, INTERP_PTS, NUM_EOS_PARAMS};
/** @brief enum NODES enumerator over tree node types */
enum NODES {
  BCO1=0, BCO2=1, BINARY, BH, NS, FIELDS, \
  STAGES, SCONTROLS, SSETTINGS, NUM_NODES
};
/** @brief enum BCO_FILEDS enumerator over BCO field types */
enum BCO_FIELDS {
  CONF=0, LAPSE, SHIFT, ENTH, LOGH, NDENS, PHI, NU, INCA, BIGA, \
  NP, KS_METRIC, KS_LAPSE, KS_K, NUM_BCO_FIELDS \
};

/** @brief enum STAGES enumerator over solver stages */
enum STAGE {
  PRE=0, FIXED_OMEGA, NOROT_BC, COROT_EQUAL, \
  TOTAL, TOTAL_BC, TOTAL_FIXED_COM, TESTING, \
  GRAV,  VEL_POT_ONLY, ECC_RED, BIN_BOOST, NUM_STAGES \
};

/** @brief enum CONTROLS enumerator over sequence controls */
enum CONTROLS {
  USE_PN, USE_FIXED_R, SEQUENCES, CHECKPOINT, MB_FIXING, \
  DELETE_SHIFT, COROT_BIN, USE_CONFIG_VARS, FIXED_GOMEGA, \
  UPDATE_INIT, USE_BOOSTED_CO, ITERATIVE_CHI, USE_FIXED_LAPSE, \
  ITERATIVE_M, RESOLVE, REGRID, SAVE_COS, CO_USE_SHELLS, NUM_CONTROLS
};

enum SEQ_SETTINGS {
  PREC, MAX_ITER, INIT_RES, FINAL_CHI, NUM_SEQ_SETTINGS
};

/**@{
  * extern definitions of maps containing the maps of strings for each parameter
  * name to the associated enumerator index
  */
extern const std::map<std::string, BIN_PARAMS> MBIN_PARAMS;
extern const std::map<std::string, BCO_PARAMS> MBCO_PARAMS;
extern const std::map<std::string, EOS_PARAMS> MEOS_PARAMS;
extern const std::map<std::string, NODES> M_REQ_NODES;
extern const std::map<std::string, NODES> MBCO;
extern const std::map<std::string, BCO_FIELDS> MBCO_FIELDS;
extern const std::map<std::string, BCO_FIELDS> MBCO_VFIELDS;
extern const std::map<std::string, BCO_FIELDS> MBCO_SFIELDS_0;
extern const std::map<std::string, BCO_FIELDS> MBCO_SFIELDS_1;
extern const std::map<std::string, STAGE> MSTAGE;
extern const std::map<std::string, CONTROLS> MCONTROLS;
extern const std::map<std::string, SEQ_SETTINGS> MSEQ_SETTINGS;
extern const std::map<std::string, STAGE> MBNSSTAGE;
extern const std::map<std::string, STAGE> MBHNSSTAGE;
extern const std::map<std::string, STAGE> MBBHSTAGE;
extern const std::map<std::string, STAGE> MBHSTAGE;
extern const std::map<std::string, STAGE> MKSBHSTAGE;
extern const std::map<std::string, STAGE> MNSSTAGE;
/**@} end extern group definition*/
/**@} end config_enums group*/
