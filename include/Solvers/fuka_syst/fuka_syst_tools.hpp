#pragma once
#include "exporter_utilities.hpp"
#include <string>
namespace FUKA_Syst_tools {

template<class space_t, class dict_t>
void export_radii(space_t & space, dict_t& vars, 
  int const dom_min, int const dom_max, std::string forestr) {
  
  int const ndom = space.get_nbr_domains();
  Scalar space_radius(space);
  space_radius.annule_hard();
  for(int d = 0; d < ndom; ++d){
    space_radius.set_domain(d) = space.get_domain(d)->get_radius();
  }
  space_radius.std_base();

  uint cnt = 1;
  for(int i = dom_min; i < dom_max; ++i) {
    Index pos(space.get_domain(i)->get_radius().get_conf().get_dimensions());
    pos.set(0) = space.get_domain(i)->get_nbr_points()(0)-1;
    vars[forestr+std::to_string(cnt)] = space_radius(i)(pos) ;
    cnt++;
  }
}

// addition of raw tensor components
template<class dict_t>
inline void dict_add_tensor_cmp(System_of_eqs& syst, dict_t& vars,
  std::string var, Tensor field) {
  int c = 0;
  vars[var.c_str()]  = syst.give_val_def(("Trace"+var).c_str());
  for(std::string coord : {"_11", "_12", "_13", "_22", "_23", "_33"}) {
    auto tidx = export_utils::R2TensorSymmetricIndices[c];
    Array<int> ind (field.indices(tidx));
    vars[(var+coord).c_str()] = field(ind);
    c++;
  }
};

template<class dict_t>
inline void dict_add_vector_cmp(System_of_eqs& syst, dict_t& vars,
  std::string var, Tensor field) {
  int c = 0;
  for(std::string coord : {"_1", "_2", "_3"}) {
    Array<int> ind (field.indices(c));
    vars[(var+coord).c_str()] = field(ind);
    c++;
  }
};

inline 
std::vector<int> vector_of_domains(int const dom_min, int const dom_max) {
  #include <numeric>
  assert(dom_min <= dom_max);
  // dom_max is not added to the list
  std::vector<int> doms(dom_max-dom_min);
  std::iota(doms.begin(), doms.end(), dom_min);
  return doms;
}

}