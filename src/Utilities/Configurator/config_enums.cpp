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
  {"com",COM},                  // Center of mass shift along X axis
  {"comy",COMY},                // COM along y axis
  {"qpig",QPIG},                // units scaling - 4*pi*G
  {"rext", REXT},               // fixed exterior radius (~2*DIST)
  {"q", Q},                     // Mass ratio
  {"adot", ADOT},               // Radial infall velocity (for eccentricity reduction) 
  {"ecc_omega", ECC_OMEGA},     // Fixed omega used for eccentricity reduction 
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
  {"velx",BVELX},               // Boost along X - Fix Px
  {"vely",BVELY},               // Boost along Y - Fix Py
  {"decay_limit", DECAY},       // Decay limit to use when importing BCOs into binary
  {"kerr_chi", KERR_CHI},       // Kerr parameter a=J/M
  {"kerr_mch", KERR_MCH},       // Mass given to the analytical kerr background
  {"n_inner_shells",NINSHELLS}, // Shells inside a NS - binary only
};

const std::map<std::string, EOS_PARAMS> MEOS_PARAMS = {
  {"eostype",EOSTYPE}, // Cold_PWPoly, Cold_Table
  {"eosfile",EOSFILE}, // Polytrope or Tabulated EOS file
  {"h_cut",HCUT}, // value to cut the specific enthalpy (0 is default)
  {"interpolation_pts", INTERP_PTS} // number of points to use for interpolating the table
};

//required independent of binary, BCO, etc
const std::map<std::string, NODES> M_REQ_NODES = {
  {"fields",FIELDS},
  {"stages",STAGES},
  {"sequence_controls",SCONTROLS},
  {"sequence_settings",SSETTINGS},
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
  {"np", NP},
  {"ks_metric", KS_METRIC},
  {"ks_lapse", KS_LAPSE},
  {"ks_k", KS_K}
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
  {"binary_boost", BIN_BOOST},
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
  {"pre",PRE},                //deprecated v1 only
  {"fixed_omega",FIXED_OMEGA},//deprecated v1 only
  {"corot_equal",COROT_EQUAL},//deprecated v1 only
  {"total",TOTAL},            //deprecated v1 only
  {"total_bc",TOTAL_BC},
  {"ecc_red", ECC_RED},
};
const std::map<std::string, STAGE> MBHNSSTAGE = {
  {"total",TOTAL},
  {"total_bc",TOTAL_BC},
  {"ecc_red", ECC_RED},
};
const std::map<std::string, STAGE> MBHSTAGE = {
  {"pre",PRE},                //deprecated v1 only
  {"total",TOTAL},            //deprecated v1 only
  {"total_bc",TOTAL_BC},
  {"binary_boost", BIN_BOOST},
};
const std::map<std::string, STAGE> MNSSTAGE = {
  {"pre",PRE},                //deprecated v1 only
  {"norot_bc",NOROT_BC},
  {"total_bc",TOTAL_BC},
  {"binary_boost", BIN_BOOST},
};

const std::map<std::string, CONTROLS> MCONTROLS = {
  {"use_pn", USE_PN},            ///< Use PN eccentricity parameters - replaces ADOT and ECC_OMEGA
  {"sequences", SEQUENCES},      ///< Enable sequence generation - placeholder
  {"checkpoint", CHECKPOINT},    ///< Disable to only output after each solver stage is successful
  {"use_fixed_r", USE_FIXED_R},  ///< Solve BCO based on FIXED_R instead of Mirr, MADM, MB, etc.
  {"fixed_mb", MB_FIXING},      ///< For an isolated NS, fix using Baryonic mass
  {"delete_shift", DELETE_SHIFT},///< at the start of the solver, choose to delete the shift
  {"corot_binary", COROT_BIN},   ///< control whether a binary is purely corotating
  
  // Control whether codes such as increase resolution make updates from the config file
  // variables or directly from the numerical space
  //{"use_config_vars", USE_CONFIG_VARS},
  {"fixed_bin_omega", FIXED_GOMEGA}, ///< Fix binary orbital frequency
  {"update_initial", UPDATE_INIT}, ///< historical: add initial section to config
  {"use_boosted_co", USE_BOOSTED_CO}, ///< use boosted compact objects to construct binary initial guess
  //{"iterative_chi", ITERATIVE_CHI},
  {"fixed_lapse", USE_FIXED_LAPSE}, ///< Use fixed lapse BC on black holes
  {"resolve", RESOLVE}, ///<Force resolve of ID even if a checkpoint exists
  {"initial_regrid", REGRID}, ///< Regrid before solving from a previous solution
  {"centralized_cos", SAVE_COS},///< Save CO solutions to a central location for reuse
  {"co_use_shells", CO_USE_SHELLS}, ///< Isolated Compact objects use defined shells (binary solvers)
};

const std::map<std::string, SEQ_SETTINGS> MSEQ_SETTINGS = {
  {"solver_precision", PREC}, ///< Threshold for a converged solution
  {"solver_max_iterations", MAX_ITER}, ///< Maximum iterations before solver is terminated
  {"initial_resolution", INIT_RES}, ///< initial resolution to solve from (default 9)
};

const std::map<std::string, STAGE> MKSBHSTAGE = {
  {"total_bc",TOTAL_BC},
  {"binary_boost", BIN_BOOST},
};