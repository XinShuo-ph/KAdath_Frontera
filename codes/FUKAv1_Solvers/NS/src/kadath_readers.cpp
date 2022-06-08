/*
 * Copyright 2021
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
#include "kadath_adapted.hpp"
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"
#include "coord_fields.hpp"
#include "python_reader.hpp"
using namespace Kadath::Margherita ;


// the space type
typedef Kadath::Space_spheric_adapted space_t;

// specialized quantities for a BNS system
struct ns_vars_t : public Kadath::vars_base_t<ns_vars_t> {};
// define the actual quantities and their order in the file!
template<> Kadath::var_vector Kadath::vars_base_t<ns_vars_t>::vars = {
                                                                        {"conf", SCALAR},
                                                                        {"lapse", SCALAR},
                                                                        {"shift", VECTOR},
                                                                        {"logh", SCALAR},
                                                                      };

class ns_reader_t : public Kadath::python_reader_t<space_t, ns_vars_t> {
  std::string config_filename;
  kadath_config_boost<BCO_NS_INFO> bconfig;

  public:
  ns_reader_t(std::string const filename) : Kadath::python_reader_t<space_t, ns_vars_t>(filename),
                                             config_filename(filename.substr(0,filename.size()-3)+"info"),
                                             bconfig(config_filename) {
    // setup eos to before calling solver
    const double h_cut = bconfig.eos<double>(HCUT);
    const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
    const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

    if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
      this->compute_defs<eos_t>();
    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
                              2000 : bconfig.eos<int>(INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      this->compute_defs<eos_t>();
    }
    else { 
      std::cerr << eos_type << " is not recognized.\n";
      std::_Exit(EXIT_FAILURE);
    }
    // end eos setup and solver
   
  }

  template <typename eos_t>
  void compute_defs() {
    Kadath::Scalar const & conf = extractField<Kadath::Scalar>("conf");
    Kadath::Scalar const & lapse = extractField<Kadath::Scalar>("lapse");
    Kadath::Vector const & shift = extractField<Kadath::Vector>("shift");
    Kadath::Scalar const & logh = extractField<Kadath::Scalar>("logh");

    Base_tensor basis(shift.get_basis());
  	Metric_flat fmet (space, basis) ;
  	int ndom = space.get_nbr_domains() ;

    double xo  = bco_utils::get_center(space, ndom-1);

  	// setup coordinate vector fields for System_of_eqs
    CoordFields<Space_spheric_adapted> cfields(space);
    vec_ary_t coord_vectors {default_co_vector_ary(space)};

    update_fields_co(cfields, coord_vectors, {}, xo);
    // end coordinate field setup

  	double loghc = bco_utils::get_boundary_val(0, logh, INNER_BC);

  	// Setup system of equations and definitions
    System_of_eqs syst (space, 0, ndom-1) ;
    fmet.set_system(syst, "f") ;

    // user defined OPEs to access EOS
    Param p;
    syst.add_ope("eps"   , &EOS<eos_t,EPSILON>::action, &p);
    syst.add_ope("press" , &EOS<eos_t,PRESSURE>::action, &p);
    syst.add_ope("rho"   , &EOS<eos_t,DENSITY>::action, &p);
    syst.add_ope("dHdlnrho"   , &EOS<eos_t,DHDRHO>::action, &p);
    
    // some constants
    syst.add_cst("4piG", bconfig(BCO_QPIG));
    syst.add_cst("PI"    , M_PI)  ;

    // coordinate fields
    syst.add_cst("mg"  , *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("ex"  , *coord_vectors[EX]);
    syst.add_cst("ey"  , *coord_vectors[EX]);
    syst.add_cst("ez"  , *coord_vectors[EX]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    // NS1 characteristics
    syst.add_cst("Mb"   , bconfig(MB)) ;
    syst.add_cst("chi"  , bconfig(CHI)) ;
    syst.add_cst("Madm" , bconfig(MADM)) ;
    syst.add_cst("ome" , bconfig(OMEGA)) ;
    syst.add_cst("Hc"   , loghc) ;

    // Fields
    syst.add_cst("P"     , conf) ;
    syst.add_cst("N"     , lapse) ;
    syst.add_cst("bet"   , shift) ;
    syst.add_cst("H"     , logh) ;

    //Useful definitions
    syst.add_def("NP = P*N");
    syst.add_def("Ntilde = N / P^6");
    
    // rotation vector field, including the shift
    syst.add_def("omega^i = bet^i + ome * mg^i");

    syst.add_def("h = exp(H)") ;
    syst.add_def("press = press(h)");
    syst.add_def("eps = eps(h)");
    syst.add_def("rho = rho(h)");
    syst.add_def("delta = h - eps - 1.");

    // Conformal Extrinsic curvature.
    syst.add_def("A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

    // ADM Linear Momentum
    syst.add_def(ndom - 1, "intPx = A_i^j * ex_j * einf^i") ;
    syst.add_def(ndom - 1, "intPy = A_i^j * ey_j * einf^i") ;
    syst.add_def(ndom - 1, "intPz = A_i^j * ez_j * einf^i") ;

    syst.add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 2. / 4piG");
    syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P * 2 / 4piG");
    syst.add_def(ndom - 1, "intMk = (einf^i * D_i N - A_ij * einf^i * bet^j) / 4piG");

      for (int d = 0; d < ndom; d++) {
        switch (d) {
        case 0:
        case 1:
          // definitions for the fluid 3-velocity
          // and its Lorentz factor
          syst.add_def(d, "U^i = omega^i / N");
          syst.add_def(d, "Usquare = P^4 * U_i * U^i");
          syst.add_def(d, "Wsquare = 1. / (1. - Usquare)");
          syst.add_def(d, "W = sqrt(Wsquare)");

          // rescaled sources and constraint equations        
          syst.add_def(d, "Etilde = press * h * Wsquare - press * delta") ;
          syst.add_def(d, "Stilde = 3 * press * delta + (Etilde + press * delta) * Usquare") ;
          syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i") ;

          syst.add_def(d, "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta + 4piG / 2. * P^5 * Etilde") ;
          syst.add_def(d, "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * A_ij *A^ij "
                                 "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
          syst.add_def(d, "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                                 "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * ptilde^i");

          // integrant of the baryonic mass integral
          syst.add_def(d, "intMb = P^6 * rho(h) * W");
          // the first integral of the Euler equation in case of a axisymmetric rotating star
          syst.add_def(d, "firstint = H + log(N) - log(W)");

          break;
        default:
          // outside the star the matter is absent and the sources are zero
          syst.add_def(d, "firstint      = h") ;

          syst.add_def(d, "eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
          syst.add_def(d, "eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
          syst.add_def(d, "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * "
                          "A^ij * D_j Ntilde");
          break;
        }
      }

    // area measurement
    syst.add_def("intMsq = P^4");

    // useful constraction definitions
    syst.add_def("drhodx = ex^i * D_i rho");
    syst.add_def("dHdx = ex^i * D_i H");
    syst.add_def("dHdx2 = ex^i * D_i dHdx");

    // add vars to boost dictionary
    vars["cP"] = syst.give_val_def("eqP");
    vars["cNP"] = syst.give_val_def("eqNP");

    Kadath::Vector beta_res(syst.give_val_def("eqbet"));

    vars["cbetx"] = beta_res(1);
    vars["cbety"] = beta_res(2);
    vars["cbetz"] = beta_res(3);

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

    Kadath::Vector v(syst.give_val_def("U"));

    vars["velx"] = v(1);
    vars["vely"] = v(2);
    vars["velz"] = v(3);

    vars["betx"] = shift(1);
    vars["bety"] = shift(2);
    vars["betz"] = shift(3);

    auto export_radii = [&](int dom, std::string forestr)
    {        
      auto radii=bco_utils::get_rmin_rmax(space, dom);
      vars[forestr+"-pole"] = radii[0] ;
      vars[forestr+"-equi"] = radii[1] ;
      
      // only valid for irrotational stars
      double A = space.get_domain(dom+1)->integ(syst.give_val_def("intMsq")()(dom+1), INNER_BC);
      vars[forestr+"-areal"]= sqrt(A / 4. / acos(-1.));
    };
    export_radii(1, "NS1_R");
  }
};

BOOST_PYTHON_MODULE(kadath_readers)
{
    // initialize python types
    Kadath::initPythonBinding<space_t>();
    Kadath::constructPythonReader<ns_reader_t>("ns_reader");
}
