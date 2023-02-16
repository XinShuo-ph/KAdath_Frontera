/*
 * Copyright 2022
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 * L. Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
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
#include "kadath_adapted_bh.hpp"
#include "bco_utilities.hpp"
#include "python_reader.hpp"
#include "include/fuka_py.hpp"

// the space type
typedef Kadath::Space_adapted_bh space_t;
using bh_config_t = kadath_config_boost<BCO_BH_INFO>;
using bh_Configurator_reader_t = Configurator_reader_t<bh_config_t>;

// specialized quantities for a BNS system
struct bh_vars_t : public Kadath::vars_base_t<bh_vars_t> {};
// define the actual quantities and their order in the file!
template<> Kadath::var_vector Kadath::vars_base_t<bh_vars_t>::vars = {
  {"conf", SCALAR},
  {"lapse", SCALAR},
  {"shift", VECTOR},
};

class bh_reader_t : public Kadath::python_reader_t<space_t, bh_vars_t> {
  std::string config_filename;
  kadath_config_boost<BCO_BH_INFO> bconfig;

  public:
  bh_reader_t(std::string const filename) 
  : Kadath::python_reader_t<space_t, bh_vars_t>(filename),
    config_filename(filename.substr(0,filename.size()-3)+"info"),
    bconfig(config_filename) { 
    
    this->compute_defs(); 
    bh_Configurator_reader_t pybconfig(config_filename);
    config = pybconfig.config;
  }

  void compute_defs() {
    Kadath::Scalar const & conf = extractField<Kadath::Scalar>("conf");
    Kadath::Scalar const & lapse = extractField<Kadath::Scalar>("lapse");
    Kadath::Vector const & shift = extractField<Kadath::Vector>("shift");

    Base_tensor basis(shift.get_basis());
  	Metric_flat fmet (space, basis) ;
  	int ndom = space.get_nbr_domains() ;

    double xo  = bco_utils::get_center(space, ndom-1);
    std::vector<int> vac_Domains{
      FUKA_Syst_tools::vector_of_domains(0, ndom)
    };

  	// setup coordinate vector fields for System_of_eqs
    CoordFields<space_t> cfields(space);
    vec_ary_t coord_vectors {default_co_vector_ary(space)};
    update_fields_co(cfields, coord_vectors, {}, xo);
    // end coordinate field setup

  	// Setup system of equations and definitions
    System_of_eqs syst (space, 0, ndom-1) ;
    fmet.set_system(syst, "f") ;

    auto init_if = [&](auto c, auto idx) {
      if(std::isnan(bconfig(idx)))
        bconfig(idx) = 0.;
      syst.add_cst(c, bconfig(idx));
    };

    // Fields - must be initialized before common setup
    syst.add_cst("P"  , conf) ;
    syst.add_cst("N"  , lapse) ;
    syst.add_cst("bet", shift) ;

    syst.add_cst("ome", bconfig(OMEGA));
    
    // Avoid excision region (d=0,1) and compactified (d=ndom-1)
    for(auto d = 2; d < ndom-1; ++d) {
      syst.add_def(d, "drP = dr(P)");
      syst.add_def(d, "ddrP = dr(drP)");
    }
    
    // Initialize constants and constant fields
    FUKA_Syst_tools::syst_init_co(syst, coord_vectors, bconfig);
    
    // Define equations in relevant domains
    FUKA_Syst_tools::syst_init_quasi_local_defs(
      syst, vac_Domains,"mg","sm"
    );
    FUKA_Syst_tools::syst_init_eqdefs_vac(syst, vac_Domains);
    FUKA_Syst_tools::syst_init_contraction_defs_vac(syst);
    
    FUKA_Syst_tools::syst_vars(vars, syst);
    FUKA_Syst_tools::syst_vars_BH(vars, syst, 2);
    FUKA_Syst_tools::export_radii(space, vars, 0, ndom-1, "BH_R");
    
    FUKA_Syst_tools::dict_add_vector_cmp(
      syst, vars, "shift", Tensor(shift)
    );
    
    // if(!std::isnan(bconfig.set(BVELX))) {
    //   syst.add_cst("xvel" , bconfig(BVELX));
    //   boosts+=" + xvel * ex^i";
    // }
    // if(!std::isnan(bconfig.set(BVELY))) {
    //   syst.add_cst("yvel" , bconfig(BVELY));
    //   boosts+=" + yvel * ey^i";
    // }
    // init_if("xboost", BVELX); 
    // init_if("yboost", BVELY);
    
    
/*    auto export_radii = [&](int dom, std::string forestr)
    {        
      auto radii=bco_utils::get_rmin_rmax(space, dom);
      vars[forestr+"-pole"] = radii[0] ;
      vars[forestr+"-equi"] = radii[1] ;
      
      // only valid for irrotational stars
      double A = space.get_domain(dom+1)->integ(syst.give_val_def("intMsq")()(dom+1), INNER_BC);
      vars[forestr+"-areal"]= sqrt(A / 4. / acos(-1.));
    };*/
    
    //create old radius scalar field
  

  }
};
BOOST_PYTHON_MODULE(_bh_reader)
{
    // initialize python types
    Kadath::initPythonBinding<space_t>();
    Kadath::constructPythonReader<bh_reader_t>("bh_reader");
}
