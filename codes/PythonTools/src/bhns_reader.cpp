/*
 * Copyright 2022
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Konrad Topolski <topolski@itp.uni-frankfurt.de>
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
#include "kadath_bhns.hpp"
#include "EOS/EOS.hh"
#include "bco_utilities.hpp"
#include "python_reader.hpp"
#include "Configurator/pyconfigurator.hpp"
#include "include/fuka_py.hpp"

using namespace Kadath::Margherita;

// the space type
typedef Kadath::Space_bhns space_t;

// specialized quantities for a BHNS system
struct bhns_vars_t : public Kadath::vars_base_t<bhns_vars_t> {};
// define the actual quantities and their order in the file!
template<> Kadath::var_vector Kadath::vars_base_t<bhns_vars_t>::vars = {
  {"conf", SCALAR},
  {"lapse", SCALAR},
  {"shift", VECTOR},
  {"logh", SCALAR},
  {"phi", SCALAR},
};

class bhns_reader_t : public Kadath::python_reader_t<space_t, bhns_vars_t> {
  std::string config_filename;
  kadath_config_boost<BIN_INFO> bconfig;

  public:

  bhns_reader_t(std::string const filename) 
    : Kadath::python_reader_t<space_t, bhns_vars_t>(filename),
      config_filename(filename.substr(0,filename.size()-3)+"info"),
      bconfig(config_filename) {
    // setup eos to before calling solver
    const double h_cut = bconfig.eos<double>(HCUT, BCO1);
    const std::string eos_file = bconfig.eos<std::string>(EOSFILE, BCO1);
    const std::string eos_type = bconfig.eos<std::string>(EOSTYPE, BCO1);

    if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
      this->compute_defs<eos_t>();
    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.eos<int>(INTERP_PTS, BCO1) == 0) ? \
                              2000 : bconfig.eos<int>(INTERP_PTS, BCO1);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      this->compute_defs<eos_t>();
    }
    else { 
      std::cerr << eos_type << " is not recognized.\n";
      std::_Exit(EXIT_FAILURE);
    }
    // end eos setup and solver
    bin_Configurator_reader_t pybconfig(config_filename);
    config = pybconfig.config;
  }

  template <typename eos_t> 
  void compute_defs() {
    Kadath::Scalar const & conf = extractField<Kadath::Scalar>("conf");
    Kadath::Scalar const & lapse = extractField<Kadath::Scalar>("lapse");
    Kadath::Vector const & shift = extractField<Kadath::Vector>("shift");
    Kadath::Scalar const & logh = extractField<Kadath::Scalar>("logh");
    Kadath::Scalar const & phi = extractField<Kadath::Scalar>("phi");

    Base_tensor basis(shift.get_basis());
    Metric_flat fmet (space, basis) ;
    int ndom = space.get_nbr_domains() ;

    double xc1 = bco_utils::get_center(space, space.NS);
    double xc2 = bco_utils::get_center(space, space.BH);
    double xo  = bco_utils::get_center(space, ndom-1);

    // setup coordinate vector fields for System_of_eqs
    CoordFields<Space_bhns> cfields(space);
    vec_ary_t coord_vectors = default_binary_vector_ary(space);
    update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);

    Vector CART(space, CON, basis);
    CART = cfields.cart();
    // end coordinate field setup

    std::vector<int> NS_int_Domains{
      FUKA_Syst_tools::vector_of_domains(space.NS, space.ADAPTEDNS+1)
    };
    std::vector<int> NS_ext_Domains{
      FUKA_Syst_tools::vector_of_domains(space.ADAPTEDNS+1, space.BH)
    };
    std::vector<int> BH_Domains{
      FUKA_Syst_tools::vector_of_domains(space.BH+2, space.OUTER)
    };
    std::vector<int> vac_Domains{
      FUKA_Syst_tools::vector_of_domains(space.ADAPTEDNS+1, ndom)
    };

    // Setup system of equations and definitions
    System_of_eqs syst (space, 0, ndom-1) ;
    fmet.set_system(syst, "f") ;

    syst.add_cst("P"  , conf) ;
    syst.add_cst("N"  , lapse) ;
    syst.add_cst("bet", shift) ;
    syst.add_cst("H"  , logh) ;
    syst.add_cst("phi", phi) ;
      
    // Binary Quantities
    syst.add_cst("ome"  , bconfig(GOMEGA));
    syst.add_cst("xaxis", bconfig(COM));
    syst.add_cst("yaxis", bconfig(COMY));

    // Component quantities
    syst.add_cst("omesm", bconfig(OMEGA, BCO1)) ;
    syst.add_cst("omesp", bconfig(OMEGA, BCO2)) ;

    for(auto d = 0; d < ndom; ++d) {
      if( d != space.BH && d != space.BH+1 && d < space.OUTER){
        syst.add_def(d, "drP = dr(P)");
        syst.add_def(d, "ddrP = dr(drP)");
      }
    }
    // Setup constants and constant fields
    FUKA_Syst_tools::syst_init_binary(syst, coord_vectors, bconfig);
    FUKA_Syst_tools::syst_init_inspiral(syst, bconfig, CART);
    
    // Define QL equations in relevant domains
    FUKA_Syst_tools::syst_init_quasi_local_defs(
      syst, NS_ext_Domains, "mm","sm"
    );
    FUKA_Syst_tools::syst_init_quasi_local_defs(
      syst, BH_Domains, "mp","sp"
    );
    
    // Setup hydro operators and definitions
    FUKA_Syst_tools::syst_init_defs_hydro<eos_t>(syst);
    // Set hydro contraction definitions
    FUKA_Syst_tools::syst_init_contraction_defs_hydro(syst);
    // Set hydro constraint definitions in relevant domains
    FUKA_Syst_tools::syst_init_eqdefs_hydro(
      syst, NS_int_Domains, "s^i  = omesm * mm^i"
    );
    // Set hydro QL definitions in relevant domains
    FUKA_Syst_tools::syst_init_quasi_local_defs_hydro(syst, NS_int_Domains);

    // Set vacuum constraints definitions in relevant domains
    FUKA_Syst_tools::syst_init_eqdefs_vac(syst, vac_Domains);
    // Set vacuum constraction definitions
    FUKA_Syst_tools::syst_init_contraction_defs_vac(syst);

    // Fill vars dictionary
    FUKA_Syst_tools::syst_vars(vars, syst);
    FUKA_Syst_tools::syst_vars_hydro(vars, syst);
    FUKA_Syst_tools::syst_vars_BH(vars, syst, space.ADAPTEDBH+1, "BH_");
    FUKA_Syst_tools::syst_vars_NS(
      vars, syst, space.ADAPTEDNS+1, bconfig(MADM, BCO1), NS_int_Domains, "NS_"
    );

    vars["xc1"] = xc1;
    vars["xc2"] = xc2;
    vars["xo"] = xo;

    FUKA_Syst_tools::dict_add_vector_cmp(syst, vars, "shift", Tensor(shift));

    FUKA_Syst_tools::export_radii(
      space, vars, space.ADAPTEDNS, space.BH, "NS_R"
    );
    FUKA_Syst_tools::export_radii(
      space, vars, space.ADAPTEDBH, space.OUTER, "BH_R"
    );
    
    Scalar space_radius(space);
    space_radius.annule_hard();
    for(int d = 0; d < ndom; ++d){
      space_radius.set_domain(d) = space.get_domain(d)->get_radius();
    }
    space_radius.std_base();
    vars["rfield"] = space_radius;

	}
};

BOOST_PYTHON_MODULE(_bhns_reader)
{
    // initialize python types
    Kadath::initPythonBinding<space_t>();
    Kadath::constructPythonReader<bhns_reader_t>("bhns_reader");
}
