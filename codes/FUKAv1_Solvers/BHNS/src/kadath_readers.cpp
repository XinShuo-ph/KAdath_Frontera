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
#include "Configurator/config_binary.hpp"
#include "coord_fields.hpp"
#include "bco_utilities.hpp"

using namespace Kadath::Margherita ;

//template <var_t var>
//using eos = EOS<Cold_PWPoly,var>;
//using eos = EOS<Cold_Table,var>;

#include "python_reader.hpp"
#include "coord_fields.hpp"

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
  bhns_reader_t(std::string const filename) : Kadath::python_reader_t<space_t, bhns_vars_t>(filename),
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
    
    double Mtot = bconfig(MADM, BCO1) + bconfig(MCH, BCO2); 

    double xc1 = bco_utils::get_center(space, space.NS);
    double xc2 = bco_utils::get_center(space, space.BH);
    double xo  = bco_utils::get_center(space, ndom-1);

    CoordFields<Space_bhns> cfields(space);
    vec_ary_t coord_vectors = default_binary_vector_ary(space);
    update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);

    Vector CART(space, CON, basis);
    CART = cfields.cart();
    // end coordinate field setup

    double loghc = bco_utils::get_boundary_val(space.NS, logh, INNER_BC);

    // Setup system of equations and definitions
    System_of_eqs syst (space, 0, ndom-1) ;
    fmet.set_system(syst, "f") ;

    // user defined OPEs to access EOS
    Param p;
    syst.add_ope ("eps"   , &EOS<eos_t,EPSILON>::action, &p);
    syst.add_ope ("press" , &EOS<eos_t,PRESSURE>::action, &p);
    syst.add_ope ("rho"   , &EOS<eos_t,DENSITY>::action, &p);
    syst.add_ope ("dHdlnrho"   , &EOS<eos_t,DHDRHO>::action, &p);
    
    // some constants
    syst.add_cst ("4piG"  , bconfig(QPIG)) ;
    syst.add_cst ("PI"    , M_PI)  ;

    // coordinate fields
    syst.add_cst ("mg"    , *coord_vectors[GLOBAL_ROT]) ;
    syst.add_cst ("mm"    , *coord_vectors[BCO1_ROT]) ;
    syst.add_cst ("mp"    , *coord_vectors[BCO2_ROT]) ;

    syst.add_cst ("ex"    , *coord_vectors[EX])  ;
    syst.add_cst ("ey"    , *coord_vectors[EY])  ;
    syst.add_cst ("ez"    , *coord_vectors[EZ])  ;

    syst.add_cst ("sm"    , *coord_vectors[S_BCO1])  ;
    syst.add_cst ("sp"    , *coord_vectors[S_BCO2])  ;
    syst.add_cst ("einf" , *coord_vectors[S_INF])  ;

    syst.add_cst ("Mb"   , bconfig(MB    , BCO1)) ;
    syst.add_cst ("chins"  , bconfig(CHI   , BCO1)) ;
    syst.add_cst ("Madmns" , bconfig(MADM, BCO1)) ;
    syst.add_cst ("omes" , bconfig(OMEGA , BCO1)) ;
    syst.add_cst ("Hc"   , loghc) ;

    // BH characteristics
    syst.add_cst    ("Mirr"   , bconfig(MIRR,BCO2)) ;
    syst.add_cst    ("chibh" , bconfig(CHI,BCO2)) ;
    syst.add_cst    ("CMbhch"  , bconfig(MCH,BCO2));

    // Binary characteristics
    syst.add_cst ("xaxis" , bconfig(COM)) ;
    syst.add_cst ("yaxis" , bconfig(COMY)) ;
    syst.add_cst ("ome"   , bconfig(GOMEGA)) ;

    std::string eccstr{};
    if(!std::isnan(bconfig.set(ADOT))) {
      syst.add_cst    ("adot" , bconfig(ADOT));
      syst.add_cst    ("r"   , CART);
      syst.add_def    ("comr^i = r^i - xaxis * ex^i + yaxis * ey^i");
      eccstr +=" + adot * comr^i";
    }

    syst.add_cst ("P"     , conf) ;
    syst.add_cst ("N"     , lapse) ;
    syst.add_cst ("bet"   , shift) ;
    syst.add_cst ("H"     , logh) ;
    syst.add_cst ("phi"   , phi) ;

    //Useful definitions
    syst.add_def ("NP = P*N");
    syst.add_def ("Ntilde = N / P^6");



   // Define COM shifted orbital rotation field
    syst.add_def ("Morb^i = mg^i + xaxis * ey^i + yaxis * ex^i");
    
    // Define global shift (inertial + corotating)
    std::string omegastr {"omega^i= bet^i + ome * Morb^i" + eccstr};
    syst.add_def (omegastr.c_str());

    syst.add_def ("h = exp(H)") ;
    syst.add_def ("press = press(h)");
    syst.add_def ("eps = eps(h)");
    syst.add_def ("rho = rho(h)");
    syst.add_def ("dHdlnrho = dHdlnrho(h)");
    syst.add_def ("delta = h - eps - 1.");

    // quasi-local spin field definitions
    for(int d = space.NS; d <= space.ADAPTEDNS; ++d){
      syst.add_def(d, "s^i  = omes * mm^i");
      syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
    }


    // Conformal Extrinsic curvature.
    syst.add_def ("A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

    // ADM Linear Momentum
    syst.add_def (ndom-1, "intPx = A_i^j * ex_j * einf^i") ;
    syst.add_def (ndom-1, "intPy = A_i^j * ey_j * einf^i") ;
    syst.add_def (ndom-1, "intPz = A_i^j * ez_j * einf^i") ;

    // Komar mass measured at inf
    syst.add_def    (ndom - 1, "Mk   =  dr(N) / 4 / PI");

    // Contractions of definitions and fields for analysis
    syst.add_def    ("Axx     = A^ij * ex_i * ex_j ");
    syst.add_def    ("Ayy     = A^ij * ey_i * ey_j");
    syst.add_def    ("Azz     = A^ij * ez_i * ez_j");
    syst.add_def    ("Axy     = A^ij * ex_i * ey_j ");
    syst.add_def    ("Ayz     = A^ij * ey_i * ez_j");
    syst.add_def    ("Axz     = A^ij * ex_i * ez_j");
    syst.add_def    ("TraceA  = A^ij * f_ij");

    syst.add_def    ("betx    = bet^i  * ex_i ");
    syst.add_def    ("bety    = bet^i  * ey_i");
    syst.add_def    ("betz    = bet^i  * ez_i");

    syst.add_def    ("omegax      = omega^i  * ex_i ");
    syst.add_def    ("omegay      = omega^i  * ey_i");
    syst.add_def    ("omegaz      = omega^i  * ez_i");


   // quasi-local spin definitions
    syst.add_def (space.ADAPTEDNS+1, "intS1 = A_ij * mm^i * sm^j / 8. / PI") ;

    for (int d=0 ; d<ndom ; d++) {
      if((d >= space.BH) || d == space.ADAPTEDNS+1){				//in vacuum the matter terms in the constraint equations are absent 
        syst.add_def(d,"eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8") ;
        syst.add_def(d,"eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
        syst.add_def(d,"eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j Ntilde");

        //set to something that is defined
        syst.add_def(d, "eqphi = h");
        syst.add_def(d, "firstint      = h") ;

      }
      else {								     //in the domains containing matter fields the constraint equations have additional terms
        // Only valid for corotation
        syst.add_def(d, "Ucor^i    = omega^i / N");
        syst.add_def(d, "Ucorsquare = P^4 * Ucor_i * Ucor^i") ;
        syst.add_def(d, "Wcorsquare = 1 / (1 - Ucorsquare)");
        syst.add_def(d, "Wcor = sqrt(Wcorsquare)");
        // end corot

        syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
        syst.add_def(d, "W      = sqrt(Wsquare)");

        syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
        syst.add_def(d, "V^i    = N * U^i - omega^i");

        syst.add_def(d, "Usquare= P^4 * U_i * U^i") ;

        // First Integral
        syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)") ;

        // velocity potential equations
        syst.add_def(d, "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i (P^6 * W * V^i)");

        // Source term defintions
        syst.add_def(d, "Etilde = press * h * Wsquare - press * delta") ;
        syst.add_def(d, "Stilde = 3 * press * delta + (Etilde + press * delta) * Usquare") ;
        syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i") ;

        // Constraint equation definitions inside the star
        syst.add_def(d, "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta + 4piG / 2. * P^5 * Etilde") ;
        syst.add_def(d, "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * A_ij *A^ij "
                               "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
        syst.add_def(d, "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                               "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * ptilde^i");

        // Baryonic Mass
        syst.add_def(d, "intMb  = P^6 * rho * W") ;
        // QLMADM
        syst.add_def(d, "intM   = - D_i D^i P * 2. / 4piG") ;
      }
    }

    // ADM Quantities
    syst.add_def    (ndom - 1, "intJ = multr(A_ij * Morb^j * einf^i) / 8 / PI");
    syst.add_def    (ndom - 1, "Madm = -dr(P) / 2 / PI");

    // Black hole spin measurement
    syst.add_def    ("intSBH   = A_ij * mp^i * sp^j   / 8 / PI") ;

    // area measurement around the black hole 
    syst.add_def ("intMsqBH  = P^4 / 16. / PI") ;

    // boost definitions for the BH
    double SBH     = space.get_domain(space.BH+2)->integ(syst.give_val_def("intSBH")()(space.BH+2) , OUTER_BC);
    double Mirr  = sqrt(space.get_domain(space.BH+2)->integ(syst.give_val_def("intMsqBH")()(space.BH+2) , INNER_BC));

    // computation
    vars["SBH"]    = SBH;
    vars["Mch"]  = std::sqrt(Mirr * Mirr + SBH * SBH / 4. / Mirr / Mirr);

   // useful constraction definitions
    syst.add_def("drhodx = ex^i * D_i rho");
    syst.add_def("dHdx = ex^i * D_i H");
    syst.add_def("dHdx2 = ex^i * D_i dHdx");

    // add vars to boost dictionary
    vars["c0"] = syst.give_val_def("eqP");
    vars["c1"] = syst.give_val_def("eqNP");

    Kadath::Vector beta_res(syst.give_val_def("eqbet"));

    vars["c2"] = beta_res(1);
    vars["c3"] = beta_res(2);
    vars["c4"] = beta_res(3);

    vars["c5"] = syst.give_val_def("eqphi");

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

    vars["Jadm"]  = space.get_domain(ndom-1)->integ(syst.give_val_def("intJ")()(ndom-1) , OUTER_BC);
    vars["Madm"]  = space.get_domain(ndom-1)->integ(syst.give_val_def("Madm")()(ndom-1) , OUTER_BC);
    vars["Mk"]    = space.get_domain(ndom-1)->integ(syst.give_val_def("Mk")()(ndom-1) , OUTER_BC);
    vars["Px"]    = space.get_domain(ndom - 1)->integ(syst.give_val_def("intPx")()(ndom-1), OUTER_BC);
    vars["Py"]    = space.get_domain(ndom - 1)->integ(syst.give_val_def("intPy")()(ndom-1), OUTER_BC);
    vars["Pz"]    = space.get_domain(ndom - 1)->integ(syst.give_val_def("intPz")()(ndom-1), OUTER_BC);

    // Conformal Aij contractions
    vars["Axx"]   = syst.give_val_def("Axx");
    vars["Ayy"]   = syst.give_val_def("Ayy");
    vars["Azz"]   = syst.give_val_def("Azz");
    vars["Axy"]   = syst.give_val_def("Axy");
    vars["Ayz"]   = syst.give_val_def("Ayz");
    vars["Axz"]   = syst.give_val_def("Axz");
    vars["A"]     = syst.give_val_def("TraceA");

    // total shift contractions
    vars["omegax"]    = syst.give_val_def("omegax");
    vars["omegay"]    = syst.give_val_def("omegay");
    vars["omegaz"]    = syst.give_val_def("omegaz");

    vars["xc1"] = xc1;
    vars["xc2"] = xc2;
    vars["xo"] = xo;

    // general area measurement
    syst.add_def ("intMsq  = P^4 / 16. / PI") ;

    auto export_radii = [&](int dom, std::string forestr)
    {        
      auto radii=bco_utils::get_rmin_rmax(space, dom);
      vars[forestr+"-pole"] = radii[0] ;
      vars[forestr+"-equi"] = radii[1] ;

      // only valid for irrotational COs
      double A = space.get_domain(dom+1)->integ(syst.give_val_def("intMsq")()(dom+1), INNER_BC);
      vars[forestr+"-areal"]= sqrt(A / 4. / acos(-1.));
    };
    export_radii(space.ADAPTEDNS, "NS_R");
    export_radii(space.ADAPTEDBH, "BH_R");
	}
};

BOOST_PYTHON_MODULE(kadath_readers)
{
    // initialize python types
    Kadath::initPythonBinding<space_t>();
    Kadath::constructPythonReader<bhns_reader_t>("bhns_reader");
}

