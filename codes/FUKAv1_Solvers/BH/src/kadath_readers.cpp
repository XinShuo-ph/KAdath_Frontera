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
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include "python_reader.hpp"

// the space type
typedef Kadath::Space_adapted_bh space_t;

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
  bh_reader_t(std::string const filename) : Kadath::python_reader_t<space_t, bh_vars_t>(filename),
                                             config_filename(filename.substr(0,filename.size()-3)+"info"),
                                             bconfig(config_filename) { this->compute_defs(); }

  void compute_defs() {
    Kadath::Scalar const & conf = extractField<Kadath::Scalar>("conf");
    Kadath::Scalar const & lapse = extractField<Kadath::Scalar>("lapse");
    Kadath::Vector const & shift = extractField<Kadath::Vector>("shift");

    Base_tensor basis(shift.get_basis());
  	Metric_flat fmet (space, basis) ;
  	int ndom = space.get_nbr_domains() ;

    double xo  = bco_utils::get_center(space, ndom-1);

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
    
    // some constants
    syst.add_cst("4piG", bconfig(BCO_QPIG));
    syst.add_cst("PI"  , M_PI)  ;

    // coordinate fields
    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EX]);
    syst.add_cst("ez", *coord_vectors[EX]);
    syst.add_cst("sm", *coord_vectors[S_BCO1]) ;

    syst.add_cst("ex"  , *coord_vectors[EX])  ;
    syst.add_cst("ey"  , *coord_vectors[EY])  ;
    syst.add_cst("ez"  , *coord_vectors[EZ])  ;
    syst.add_cst("einf", *coord_vectors[S_INF]) ;

    // BH characteristics
    syst.add_cst("M"    , bconfig(MIRR)) ;
    syst.add_cst("chi"  , bconfig(CHI)) ;
    syst.add_cst("CM"   , bconfig(MCH));
    syst.add_cst("ome"  , bconfig(OMEGA));
    init_if("xboost", BVELX); 
    init_if("yboost", BVELY);

    // Fields
    syst.add_cst("P"     , conf) ;
    syst.add_cst("N"     , lapse) ;
    syst.add_cst("bet"   , shift) ;

    //Useful definitions
    syst.add_def("NP = P*N");
    syst.add_def("Ntilde = N / P^6");
    
    // rotation vector field, including the shift
    syst.add_def("omega^i = bet^i + ome * mg^i");

    // Conformal Extrinsic curvature.
    syst.add_def("A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

    // ADM Linear Momentum
    syst.add_def(ndom - 1, "intPx = A_i^j * ex_j * einf^i / 8 / PI") ;
    syst.add_def(ndom - 1, "intPy = A_i^j * ey_j * einf^i / 8 / PI") ;
    syst.add_def(ndom - 1, "intPz = A_i^j * ez_j * einf^i / 8 / PI") ;

    syst.add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 2. / 4piG");
    syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P * 2 / 4piG");
    syst.add_def(ndom - 1, "intMk = (einf^i * D_i N - A_ij * einf^i * bet^j) / 4piG");
    
    syst.add_def("intS = A_ij * mg^i * sm^j / 8. / PI") ;
    syst.add_def("intMsq = P^4 / 16. / PI") ;

    syst.add_def("eqP = D^i D_i P + A_ij * A^ij / P^7 / 8") ;
    syst.add_def("eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def("eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j Ntilde");
    for(auto d = 0; d < ndom; ++d)
      //if(d < ndom-1)
        syst.add_def(d, "Pratio = P / dr(P)");
      //else
        //syst.add_def(d,"Pratio = P");

    // add vars to boost dictionary
    vars["cP"] = syst.give_val_def("eqP");
    vars["cNP"] = syst.give_val_def("eqNP");
    vars["Pratio"] = syst.give_val_def("Pratio");

    Kadath::Vector beta_res(syst.give_val_def("eqbet"));

    vars["cbetx"] = beta_res(1);
    vars["cbety"] = beta_res(2);
    vars["cbetz"] = beta_res(3);

    vars["betx"] = shift(1);
    vars["bety"] = shift(2);
    vars["betz"] = shift(3);

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
    Scalar space_radius(space);
    space_radius.annule_hard();
    for(int d = 0; d < ndom; ++d){
      space_radius.set_domain(d) = space.get_domain(d)->get_radius();
    }
    space_radius.std_base();
    auto export_radii = [&](int dom_min, int dom_max, std::string forestr)
    {
      int cnt = 1;
      for(int i = dom_min; i < dom_max; ++i) {
        Index pos(space.get_domain(i)->get_radius().get_conf().get_dimensions());
        pos.set(0) = space.get_domain(i)->get_nbr_points()(0)-1;
        vars[forestr+std::to_string(cnt)] = space_radius(i)(pos) ;
        cnt++;
      }
    };    
    export_radii(0, ndom, "BH_R");
    
    auto add_tensor_cmp = [&](std::string var, auto field) {
      int c = 0;
      for(std::string coord : {"XX", "XY", "XZ", "YY", "YZ", "ZZ"}) {
        Array<int> ind (field.indices(c));
        vars[(var+coord).c_str()] = field(ind);
        c++;
      }
    };
    auto add_vector_cmp = [&](std::string var, auto field) {
      int c = 0;
      for(std::string coord : {"X", "Y", "Z"}) {
        Array<int> ind (field.indices(c));
        vars[(var+coord).c_str()] = field(ind);
        c++;
      }
    };
    
    auto A = syst.give_val_def("A");
    add_tensor_cmp("A", A);
    add_vector_cmp("bet", Tensor(shift));
    //export_radii(1, "BH_R");
  
  }
};

BOOST_PYTHON_MODULE(kadath_readers)
{
    // initialize python types
    Kadath::initPythonBinding<space_t>();
    Kadath::constructPythonReader<bh_reader_t>("bh_reader");
}
