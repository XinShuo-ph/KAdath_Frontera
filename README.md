# Frankfurt University/Kadath Initial Data branch
### Author(s)    : L. Jens Papenfort, Samuel D. Tootle, Philippe Grandclément

## Overview - Updated
Included are the Frankfurt initial data solvers and utilities based on the Kadath
  spectral solver library.  The original solvers written by the aformentioned authors
  (hereafter denoted as FUKAv1) are located in ./codes/FUKAv1/[BH, NS, BHNS, BNS, BBH] respectively.
  The solvers in the next version release, denoted as FUKAv2, can be found in ./codes/FUKAv2/[BH, NS, BHNS, BNS, BBH] respectively.
	Both FUKAv1 (v1) and FUKAv2 (v2) includes support for polytropic equations of state as well as tabulated EOS in
  in the standard LORENE format.  Examples and additional details can be found in the [eos](./eos/) directory.

## FUKAv2 Notes: 
The release of FUKAv2 is considerable step forward in reliable generation of extremal spinning, asymmetric
binary initial data using the KADATH spectral software.  The v2 solvers aim to maximize convergence by using super-imposed
isolated solutions to setup the initial guess for binary ID.  Additionally, v2 aims to automate the generation of ID by minimizing
the workflow for the user.  Finally, v2 allows quite for considerable flexibility in setting up the config file to make ID generation 
as efficient as possible within the KADATH framework including reusing previously solved implicit isolated solutions.  For more details,
please see the documentation in the [FUKAv2](./codes/FUKAv2/).

  - FUKAv1 Specific: There have been some core changes to various utilities used when constructing ID that have required some refactoring of the v1 solvers.  Spot testing has been done to ensure these codes function as originally intended - however - these codes, overall, remain unchanged.
  
## Maintainer(s):  

Samuel D. Tootle - tootle@itp.uni-frankfurt.de (primary),  

L. Jens Papenfort - papenfort@itp.uni-frankfurt.de  

License      : GPLv3+ for all other code  

# 1. Purpose

This collection of ID solvers aims at delivering consistent initial data (ID)  
solutions to the eXtended Conformal Thin-Sandwich (XCTS) formulation of  
Einstein's field equations for a variety of compact object configurations. 
  
As each solver has their own specific nuances and considerations, we have included  
a README in each solver directory to provide a basis for getting started with the  
respective solver.  
  
Additionally, each initial data has a respective exporter which can be seen  
in src/Utilities/Exporters.  These exporters allow one to compile an interface code  
for an evolution toolkit along with the Kadath static library located in ./lib, in  
order to export data based on input grid points.  

# 2. Modifications from base Kadath

In addition to the solving routines included within ./codes, we also note the major overall modifications  
and additions that differ from base Kadath.  
1.  This branch includes memory optimizations that inspired portions of the optimization branch  
2.  Modification/addition of numerical spaces for the BH, BBH, BNS, and BHNS  
3.  Addition of an equation of state infrastructure utilizing Margherita standalone to handle
tabulated and polytropic EOS - see [include/EOS](./include/EOS/)  
4.  Addition of the Configurator framework to enable extensibility of solvers by managing controls,
stages, and key variables - see [include/Configurator](./include/Configurator)  
5.  Addition of exporters for all the previously mentioned ID types - see [src/Utilities/Exporters](./include/Configurator)

Note: as of summer 2021, the FUKA solvers are based on the deprecated branch of Kadath.  Given the optimizations and changes made
within the FUKA branch conflict with those implimented in the `master` branch (previously the `optimized` branch), a considerable
level of effort is required to merge FUKA with the new `master` branch as well as test to see which optimizations provide better results.
Currently, there is no timeline for when this will be done.

# 3. Public Thorns for use with the Einstein Toolkit

The following workspace includes the FUKA ID respository (including versioned branches) 
as well as available thorns for use with the Einstein Toolkit in order to import FUKA ID:
https://bitbucket.org/fukaws/

# REQUIRED CITATIONS:

1) L. Jens Papenfort, Samuel D. Tootle, Philippe Grandclément, Elias R. Most, Luciano Rezzolla: https://arxiv.org/abs/2103.09911  
  
2) Philippe Grandclément, http://dx.doi.org/10.1016/j.jcp.2010.01.005  
