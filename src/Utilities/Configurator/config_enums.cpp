/*
    Copyright 2019 Samuel Tootle

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

#include <Configurator/config_enums.hpp>
#include <string>
#include <map>

const std::map<std::string, BIN_PARAMS> MBIN_PARAMS = {
  {"res",BIN_RES},              // Resolution
  {"distance",DIST},            // Separation distance
  {"d_dist",DDIST},             // Change in separation distance
  {"global_omega",GOMEGA},      // Orbital angular velocity
  {"velx",BVELX},               // Boost along X - Fix Px
  {"vely",BVELY},               // Boost along Y - Fix Py
  {"com",COM},                  // Center of mass shift along X axis
  {"comy",COMY},                // COM along y axis
  {"qpig",QPIG},                // units scaling - 4*pi*G
  {"rext", REXT},               // fixed exterior radius (~2*DIST)
  {"q", Q},                     // Mass ratio
  {"adot", ADOT},                
  {"ecc_omega", ECC_OMEGA},     // Fixed omega used for eccentricity reduction                
  {"init_res",INIT_RES},        // Initial Resolution to use before res increase
  {"outer_shells", OUTER_SHELLS}, // Number of shells before compactified domain
};

const std::map<std::string, BCO_PARAMS> MBCO_PARAMS = {
  {"res",BCO_RES},              // Resolution
  {"qpig",BCO_QPIG},            // units scaling - 4*pi*G
  {"rin",RIN},                  // nucleus fixed radius
  {"rmid",RMID},                // surface radius guess
  {"fixed_r",FIXED_R},          // fixed surface radius
  {"rout",ROUT},                // fixed outer domain radius
  {"nshells",NSHELLS},
  {"omega",OMEGA},              // Angular velocity
  {"chi",CHI},                  // Dimensionless spin
  {"mirr",MIRR},                // Irreducible Mass
  {"mch",MCH},                  // Christodoulou Mass
  {"mb",MB},                    // Baryonic Mass
  {"nc",NC},                    // Central Density
  {"hc",HC},                    // Central Enthalpy
  {"fixed_lapse",FIXED_LAPSE},  // fixed lapse on the BH horizon
  {"madm",MADM},                // MADM of the isolated NS
  {"ql_madm",QLMADM},           // quasi-local MADM from the BNS solver
  {"dim",DIM},
  {"use_tov1d", USE_TOV1D},      // Use 1D TOV estimates for R, NC, and HC
  {"fixed_omega", FIXED_BCOMEGA},// fixed Angular velocity - chi is ignored
};

const std::map<std::string, EOS_PARAMS> MEOS_PARAMS = {
  {"eostype",EOSTYPE},
  {"eosfile",EOSFILE},
  {"h_cut",HCUT},
  {"interpolation_pts", INTERP_PTS}
};

//required independent of binary, BCO, etc
const std::map<std::string, NODES> M_REQ_NODES = {
  {"fields",FIELDS},
  {"stages",STAGES},
  {"sequence_controls",SCONTROLS}
};

const std::map<std::string, NODES> MBCO = {
  {"bh",BH},
  {"ns",NS}
};

// Fields used in Configurator
const std::map<std::string, BCO_FIELDS> MBCO_FIELDS = {
  {"conf", CONF},
  {"lapse", LAPSE},
  {"shift", SHIFT},
  {"enth", ENTH},
  {"logh", LOGH},
  {"ndens", NDENS},
  {"phi", PHI},
  {"nu", NU},
  {"incA", INCA},
  {"bigA", BIGA},
  {"np", NP}
};

// Subset of fields that are Scalars and are initialized to 0
const std::map<std::string, BCO_FIELDS> MBCO_SFIELDS_0 = {
  {"logh", LOGH},
  {"ndens", NDENS},
};

// Subset of fields that are Scalars and are initialized to 1
const std::map<std::string, BCO_FIELDS> MBCO_SFIELDS_1 = {
  {"conf", CONF},
  {"lapse", LAPSE},
  {"enth", ENTH},
  {"nu", NU},
  {"incA", INCA},
  {"bigA", BIGA},
  {"np", NP}
};

// Subset of fields that are Vectors - initialized to 0 always
const std::map<std::string, BCO_FIELDS> MBCO_VFIELDS = {
  {"shift", SHIFT},
  {"phi", PHI},
};

// all reserved stage names
const std::map<std::string, STAGE> MSTAGE = {
  {"pre",PRE},
  {"norot_bc",NOROT_BC},
  {"fixed_omega",FIXED_OMEGA},
  {"corot_equal",COROT_EQUAL},
  {"total",TOTAL},
  {"total_bc",TOTAL_BC},
  {"total_fixed_com",TOTAL_FIXED_COM},
  {"grav",GRAV},
  {"vel_pot_only",VEL_POT_ONLY},
  {"ecc_red", ECC_RED},
  {"testing",TESTING}
};

/**
 * The following stage maps are sub-sets of MSTAGE.
 * This is used to allow only the relevant stage names for a given
 * solver to be shown in the configurator file.  However, 
 * for new solvers that are not composed of NSs or BHs and have
 * not been assigned such a subset in config_bin.hpp,
 * the default is the full list of stages.
 */
const std::map<std::string, STAGE> MBNSSTAGE = {
  {"total",TOTAL},
  {"total_bc",TOTAL_BC},
  {"ecc_red", ECC_RED},
};
const std::map<std::string, STAGE> MBBHSTAGE = {
  {"pre",PRE},
  {"fixed_omega",FIXED_OMEGA},
  {"corot_equal",COROT_EQUAL},
  {"total",TOTAL},
  {"total_bc",TOTAL_BC},
  {"ecc_red", ECC_RED},
};
const std::map<std::string, STAGE> MBHNSSTAGE = {
  {"total",TOTAL},
  {"total_bc",TOTAL_BC},
  {"ecc_red", ECC_RED},
};
const std::map<std::string, STAGE> MBHSTAGE = {
  {"pre",PRE},
  {"total",TOTAL},
  {"total_bc",TOTAL_BC}
};
const std::map<std::string, STAGE> MNSSTAGE = {
  {"pre",PRE},
  {"norot_bc",NOROT_BC},
  {"total_bc",TOTAL_BC}
};

const std::map<std::string, CONTROLS> MCONTROLS = {
  {"use_pn", USE_PN},            ///< Use PN eccentricity parameters - replaces ADOT and ECC_OMEGA
  {"sequences", SEQUENCES},      ///< Enable sequence generation - placeholder
  {"checkpoint", CHECKPOINT},    ///< Disable to only output after each solver stage is successful
  {"use_fixed_r", USE_FIXED_R},  ///< Solve BCO based on FIXED_R instead of Mirr, MADM, MB, etc.
  {"mb_fixing", MB_FIXING},      ///< For an isolated NS, fix using Baryonic mass
  {"delete_shift", DELETE_SHIFT},///< at the start of the solver, choose to delete the shift
  {"corot_binary", COROT_BIN},   ///< control whether a binary is purely corotating
  
  // Control whether codes such as increase resolution make updates from the config file
  // variables or directly from the numerical space
  {"use_config_vars", USE_CONFIG_VARS},
  {"fixed_bin_omega", FIXED_GOMEGA},
  {"update_initial", UPDATE_INIT},
};
