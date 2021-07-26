# Frankfurt University/Kadath Initial Data branch
## Author(s)    : L. Jens Papenfort, Samuel D. Tootle, Philippe Grandclément
Note : Frankfurt initial data solvers and utilities based on the Kadath  
  spectral solver library.  The included solvers written by the aformentioned authors  
  are located in ./codes/[BH, NS, BHNS, BNS, BBH] respectively.  The current release  
	includes support for polytropic equations of state as well as tabulated EOS in  
  in the standard LORENE format.  Examples can be found in the ./eos directory  
  
Maintainer(s):  
Samuel D. Tootle - tootle@itp.uni-frankfurt.de (primary),  
L. Jens Papenfort - papenfort@itp.uni-frankfurt.de  
License      : GPLv3+ for all other code  
--------------------------------------------------------------------------

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

In addition to the solving routines included ./codes, we also note the major overall modifications  
and additions that differ from base Kadath.  
1.  This branch includes memory optimizations that inspired portions of the optimization branch  
2.  Modification/addition of numerical spaces for the BH, BBH, BNS, and BHNS  
3.  Addition of an equation of state infrastructure utilizing Margherita standalone to handle
tabulated and polytropic EOS - see include/EOS  
4.  Addition of the Configurator framework to enable extensibility of solvers by managing controls,
stages, and key variables - see include/Configurator  
5.  Addition of exporters for all the previously mentioned ID types - see src/Utilities/Exporters  

# 3. Public Thorns for use with the Einstein Toolkit

The following workspace includes the publicly available thorns for use with the Einstein Toolkit
in order to import FUKA ID:
https://bitbucket.org/fukaws/

# REQUIRED CITATIONS:

1) L. Jens Papenfort, Samuel D. Tootle, Philippe Grandclément, Elias R. Most, Luciano Rezzolla: https://arxiv.org/abs/2103.09911  
  
2) Philippe Grandclément, http://dx.doi.org/10.1016/j.jcp.2010.01.005  
