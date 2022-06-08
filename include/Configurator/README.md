# Configurator

## Author(s) : Samuel D. Tootle

Maintainer(s):  
Samuel D. Tootle - tootle@itp.uni-frankfurt.de
License      : GPLv3+ for all other code 
--------------------------------------------------------------------------

# 1. Purpose

The Configurator (config) code was originally created to manage save a few variable scalars (`omega`) and fixing parameters (e.g. `mch`, `chi`), however, it has evolved to being an integral interface code to manage the construction of initial data.

The config code is based on [Boost](https://www.boost.org/) `ptree` library which allows for the creation, storage, and retrieval of a "tree" structure.  Currently, the used format is the proprietary `info` format which allows the use of comments using a semi-colon `;`.  With that said, they also support the use of XML and JSON, therefore, this change is in principle possible within the current framework.

However, storage of data requires some form of mapping in order to store and retrieve useful data quickly and in a meaningful way.  To this end the `$HOME_KADATH/include/Configurator/config_enums.hpp` repository was created to contain enums of physically motivated indicies which can be used for accessing array elements storing the related data.  But this information must also be written to file which means a meaningful map between an enum and string is required.  Such maps are then stored in `$HOME_KADATH/src/Utilities/Configurator/config_enums.cpp`.  This combination allows for easily extending the enum and map dictionaries very easily and without needing to refactor old code while, at the same time, providing physically motivated storage access making code more readable.  For example, accessing the spin of a neutron star is as simple as `bconfig(CHI)` once the config file has been ingested.

Finally, aside from storing scalar quantities and keeping tabs on book keeping elements, the config also includes various controls (on/off switches) and sequence settings (important values) that allow considerable control over how ID is generated and is the basis for how extreme binary ID has become automated.

# 2. Parameter Details

## Binary Parameters

Here we list the relevant binary parameters that can be tuned by the user within the config file

```
  {"res",BIN_RES},              ///< Resolution
  {"distance",DIST},            ///< Separation distance
  {"global_omega",GOMEGA},      ///< Orbital angular frequency
  {"com",COM},                  ///< Center of mass shift along X axis
  {"comy",COMY},                ///< COM along y axis
  {"rext", REXT},               ///< fixed exterior radius (~2*DIST)
  {"q", Q},                     ///< Mass ratio
  {"adot", ADOT},               ///< Radial infall velocity (for eccentricity reduction) 
  {"ecc_omega", ECC_OMEGA},     ///< Fixed omega used for eccentricity reduction 
  {"outer_shells", OUTER_SHELLS}, ///< Number of shells before compactified domain
```

## Compact Object Parameters

Here we list the relevant compact object parameters

```
  {"res",BCO_RES},              ///< Resolution
  {"nshells",NSHELLS},          ///< Number of spherical shells that can be added around the compact object for additional local resolution
  {"omega",OMEGA},              ///< Angular frequency
  {"chi",CHI},                  ///< Dimensionless spin
  {"mirr",MIRR},                ///< Irreducible Mass
  {"mch",MCH},                  ///< Christodoulou Mass
  {"mb",MB},                    ///< Baryonic Mass
  {"nc",NC},                    ///< Central Density
  {"hc",HC},                    ///< Central specific Enthalpy
  {"fixed_lapse",FIXED_LAPSE},  ///< fixed lapse on the BH horizon
  {"madm",MADM},                ///< MADM of the isolated NS
  {"ql_madm",QLMADM},           ///< quasi-local MADM from the BNS solver
  {"dim",DIM},                  ///< Dimension of the grid (2D/3D)
  {"fixed_omega", FIXED_BCOMEGA},///< fixed Angular frequency - chi is ignored
  {"velx",BVELX},               ///< Boost along X - Fix Px
  {"vely",BVELY},               ///< Boost along Y - Fix Py
  {"decay_limit", DECAY},       ///< Decay limit to use when importing BCOs into binary
  {"n_inner_shells",NINSHELLS}, ///< Shells inside a NS - binary only
```

## EOS Parameters

```
  {"eostype",EOSTYPE}, ///< Options are currently: Cold_PWPoly, Cold_Table
  {"eosfile",EOSFILE}, ///< Polytrope or Tabulated EOS file
  {"h_cut",HCUT}, ///< value to cut the specific enthalpy (0 is default)
  {"interpolation_pts", INTERP_PTS} ///< number of points to use for interpolating the table
```

## Controls

```
  {"use_pn", USE_PN}, ///< Use 3.5PN eccentricity parameters - replaces ADOT and ECC_OMEGA
  {"sequences", SEQUENCES},      ///< Enable generation of binary ID from scratch
  {"checkpoint", CHECKPOINT},    ///< Disable to only output after each solver stage is successful
  {"fixed_mb", MB_FIXING},       ///< For an isolated NS, fix using Baryonic mass
  {"delete_shift", DELETE_SHIFT},///< at the start of the solver, choose to delete the shift
  {"corot_binary", COROT_BIN},   ///< control whether a binary is purely corotating
  {"fixed_bin_omega", FIXED_GOMEGA},   ///< Fix binary orbital frequency
  {"update_initial", UPDATE_INIT},     ///< historical: add initial section to config
  {"use_ boosted_co", USE_BOOSTED_CO}, ///< use boosted compact objects to construct binary initial guess
  {"fixed_lapse", USE_FIXED_LAPSE},    ///< Use fixed lapse BC on black holes
  {"resolve", RESOLVE},       ///<Force resolve of ID even if a checkpoint exists
  {"initial_regrid", REGRID}, ///< Regrid before solving from a previous solution - placeholder
  {"centralized_cos", SAVE_COS},///< Save compact object solutions to a central location for reuse
  {"co_use_shells", CO_USE_SHELLS}, ///< Isolated Compact objects use defined shells (binary solvers)
```