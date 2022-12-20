#pragma once
#include "Solvers/fuka_syst/fuka_syst_tools.hpp"
#include "bco_utilities.hpp"
namespace FUKA_Syst_tools {

template<class dict_t>
void syst_vars_hydro(dict_t& vars, System_of_eqs & syst) {
  #ifdef DEBUG
    std::cout << "Loading hydro vars into dictionary.\n";
  #endif
  {
    char n[] = "eqphi";
    char* name;
    name = new char[LMAX];
    trim_spaces(name, n);
    int which = -1 ;
    int valence;
    char* name_ind = 0x0  ;
    Array<int>* type_ind = 0x0 ;
    bool found = false ;

    if(syst.isdef(0, name, which, valence, name_ind, type_ind))
      vars["cPhi"] = syst.give_val_def("eqphi");
  }
  
  vars["rho"] = syst.give_val_def("rho");
  vars["eps"] = syst.give_val_def("eps");
  vars["press"] = syst.give_val_def("press");
  vars["firstint"] = syst.give_val_def("firstint");

  vars["drhodx"] = syst.give_val_def("drhodx");
  vars["dHdx"] = syst.give_val_def("dHdx");
  vars["delta"] = syst.give_val_def("delta");

  vars["W"] = syst.give_val_def("W");
  vars["h"] = syst.give_val_def("h");
  vars["Etilde"] = syst.give_val_def("Etilde");
  vars["Stilde"] = syst.give_val_def("Stilde");
  FUKA_Syst_tools::dict_add_vector_cmp(
    syst, vars, "vel", syst.give_val_def("U")
  );
}

/**
 * @brief Fill vars dict with Neutron star details
 * 
 * @tparam dict_t boost dictionary type
 * @param vars dictionary
 * @param syst system of equations
 * @param adapt_d INNER adapted domain index
 * @param Madm ADM mass to use when computing chi
 * @param iden string identifier for the dictionary
 */
template<class dict_t>
void syst_vars_NS(dict_t& vars, System_of_eqs & syst, 
  int const adapt_d, double const Madm, std::vector<int> doms,
    std::string iden = {}) {
    #ifdef DEBUG
      std::cout << "Loading NS vars into dictionary.\n";
    #endif
    auto& space = syst.get_space();

    std::vector<double> baryonic_mass{};
    std::vector<double> int_H{};
    std::vector<double> qladm_mass{};
    for(auto& d : doms){
      baryonic_mass.push_back(syst.give_val_def("intMb")()(d).integ_volume());
      int_H.push_back(syst.give_val_def("intH")()(d).integ_volume());
      qladm_mass.push_back(syst.give_val_def("intqlMadm")()(d).integ_volume());
    }
    vars[iden+"Mb"] = std::accumulate(baryonic_mass.begin(),baryonic_mass.end(),0.);
    vars[iden+"qlMADM"] = std::accumulate(qladm_mass.begin(),qladm_mass.end(),0.);
    vars[iden+"Hvol"] = std::accumulate(int_H.begin(), int_H.end(),0.);

    auto radii=bco_utils::get_rmin_rmax(space, adapt_d-1);
    vars[iden+"coordR-pole"] = radii[0] ;
    vars[iden+"coordR-equi"] = radii[1] ;
    
    // Local angular momentum on the AH
    double S = space.get_domain(adapt_d)->integ(
      syst.give_val_def("intS")()(adapt_d), OUTER_BC
    );
    vars[std::string{iden+"S"}.c_str()] = S;
    
    // Proper area A = 4 pi r^2 = integ(intArea)
    double Area = space.get_domain(adapt_d)->integ(
        syst.give_val_def("intArea")()(adapt_d) , INNER_BC
      );
    double R = sqrt(Area / 4. / M_PI);

    // Area radius
    vars[std::string{iden+"ArealR"}.c_str()] = R;

    double Chi = S / Madm / Madm;
    vars[std::string{iden+"Chi"}.c_str()] = Chi;
}
}
