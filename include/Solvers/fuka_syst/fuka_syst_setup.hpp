#pragma once
#include "coord_fields.hpp"
#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
namespace FUKA_Syst_tools {
/**
 * \addtogroup Syst_tools
 * \ingroup FUKA
 * Centralized tools to populate and manipulate System_of_eqs
 * @{*/

template<class cfields_t>
void syst_init_csts(System_of_eqs & syst, cfields_t& coord_vectors) {
  #ifdef DEBUG
    std::cout << "Loading constants and constant fields.\n";
  #endif
  syst.add_cst("PI", M_PI);
  
  // coordinate fields
  syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);

  syst.add_cst("ex"  , *coord_vectors[EX])  ;
  syst.add_cst("ey"  , *coord_vectors[EY])  ;
  syst.add_cst("ez"  , *coord_vectors[EZ])  ;
  syst.add_cst("einf", *coord_vectors[S_INF]) ;
}

inline
void syst_init_defs(System_of_eqs & syst) {
  int const ndom = syst.get_space().get_nbr_domains();
  #ifdef DEBUG
    std::cout << "Loading global definitions.\n";
  #endif
  //Useful definitions
  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");
  // Conformal Extrinsic curvature.
  syst.add_def("A^ij  = (D^i bet^j + D^j bet^i \
                      - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

  // ADM Linear Momentum
  syst.add_def(ndom - 1, "intPx = A_i^j * ex_j * einf^i / 8 / PI") ;
  syst.add_def(ndom - 1, "intPy = A_i^j * ey_j * einf^i / 8 / PI") ;
  syst.add_def(ndom - 1, "intPz = A_i^j * ez_j * einf^i / 8 / PI") ;

  syst.add_def(ndom - 1, "intMadm = -dr(P) / 2 / PI");
  syst.add_def(ndom - 1, "intMk =  dr(N) / 4 / PI");

  // Irreducible mass integrand
  syst.add_def("intArea = P^4 / 4piG") ;
  syst.add_def("intMirrsq  = intArea / 4") ;
}

inline
void syst_init_contraction_defs_vac(System_of_eqs & syst) {
  #ifdef DEBUG
    std::cout << "Loading vacuum contraction definitions.\n";
  #endif
  int const ndom = syst.get_space().get_nbr_domains();

  // Contractions of definitions and fields for analysis
  syst.add_def("Axx     = A^ij * ex_i * ex_j ");
  syst.add_def("Ayy     = A^ij * ey_i * ey_j");
  syst.add_def("Azz     = A^ij * ez_i * ez_j");
  syst.add_def("Axy     = A^ij * ex_i * ey_j ");
  syst.add_def("Ayz     = A^ij * ey_i * ez_j");
  syst.add_def("Axz     = A^ij * ex_i * ez_j");
  syst.add_def("TraceA  = A^ij * f_ij");

  // Contractions of definitions and fields for analysis
  syst.add_def("Bx = B^i  * ex_i ");
  syst.add_def("By = B^i  * ey_i");
  syst.add_def("Bz = B^i  * ez_i");

  syst.add_def("betx = bet^i  * ex_i ");
  syst.add_def("bety = bet^i  * ey_i");
  syst.add_def("betz = bet^i  * ez_i");
}

template<class config_t>
void syst_init_inspiral(System_of_eqs & syst, config_t& bconfig, Vector& CART) {
  int const ndom = syst.get_space().get_nbr_domains();
  std::string eccstr{};
  if(!std::isnan(bconfig.set(ADOT))) {
    syst.add_cst("adot", bconfig(ADOT));
    syst.add_cst("r", CART);
    syst.add_def("comr^i = r^i - xaxis * ex^i + yaxis * ey^i");
    eccstr +=" + adot * comr^i";
  } 

  // Define COM shifted orbital rotation field
  syst.add_def("Morb^i = mg^i + xaxis * ey^i + yaxis * ex^i");
  
  // Define global shift (inertial + corotating)
  std::string Bstr {"B^i= bet^i + ome * Morb^i" + eccstr};
  syst.add_def(Bstr.c_str());

  // ADM Angular Momentum Def
  syst.add_def(ndom - 1, "intJ = multr(A_ij * Morb^j * einf^i) / 2 / 4piG");
}

template<class cfields_t, class config_t>
void syst_init_binary(System_of_eqs & syst, 
  cfields_t& coord_vectors, config_t& bconfig) {
  #ifdef DEBUG
    std::cout << "Initializing standard binary fields, constants, and definitions.\n";
  #endif
  int const ndom = syst.get_space().get_nbr_domains();
  syst.add_cst("4piG", bconfig(QPIG)) ;
  syst_init_csts(syst, coord_vectors);
  syst_init_defs(syst);

  // remaining coordinate fields
  syst.add_cst("mm", *coord_vectors[BCO1_ROT]) ;
  syst.add_cst("mp", *coord_vectors[BCO2_ROT]) ;

  syst.add_cst("sm", *coord_vectors[S_BCO1])  ;
  syst.add_cst("sp", *coord_vectors[S_BCO2])  ;
}

inline
void syst_init_quasi_local_defs(System_of_eqs& syst,
  std::vector<int> doms, std::string rotdef, std::string surfdef) {
  #ifdef DEBUG
    std::cout << "Initializing QL spin definition.\n";
  #endif
  std::string qlspin_def = "intS = A_ij * " + rotdef + "^i * " + surfdef + "^j";
  qlspin_def += " / 2 / 4piG";
  for(auto& d : doms){
    syst.add_def(d, qlspin_def.c_str()) ;
  }
}

// FIXME boosts not considered
template<class cfields_t, class config_t>
void syst_init_co(System_of_eqs & syst, 
  cfields_t& coord_vectors, config_t& bconfig) {
  #ifdef DEBUG
    std::cout << "Initializing standard CO fields, constants, and definitions.\n";
  #endif
  int const ndom = syst.get_space().get_nbr_domains();
  syst.add_cst("4piG", bconfig(BCO_QPIG));
  syst_init_csts(syst, coord_vectors);

  syst.add_cst("sm", *coord_vectors[S_BCO1])  ;

  syst.add_def("B^i = bet^i + ome * mg^i");
  syst_init_defs(syst);

  // ADM Angular Momentum Def
  syst.add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 2. / 4piG");

  // Local spin definition
  syst.add_def("intS = A_ij * mg^i * sm^j / 2 / 4piG") ;
}

inline
void syst_init_eqdefs_vac(System_of_eqs& syst, std::vector<int> doms) {
  #ifdef DEBUG
    std::cout << "Initializing vacuum constraint equation definitions.\n";
  #endif
  for(auto& d : doms) {
    syst.add_def(d,"eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8") ;
    syst.add_def(d,"eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(d,"eqbet^i = D_j D^j bet^i \
                          + D^i D_j bet^j / 3. - 2. * A^ij * D_j Ntilde");
  }
}
/** @}*/
}
#include "fuka_syst_setup_hydro.hpp"

// Constraint equations
