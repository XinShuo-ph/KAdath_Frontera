#include <kadath_bin_ns.hpp>
#include <EOS/EOS.hh>
#include <coord_fields.hpp>
#include <Configurator/config_binary.hpp>
#include <exporter_utilities.hpp>
#include <bco_utilities.hpp>

#include <cmath>

using namespace Kadath;
using namespace export_utils;

std::array<std::vector<double>,NUM_OUT> KadathExportNS(int const npoints,
                                                        double const * xx, double const * yy, double const * zz,
                                                        char const * fn) {
  kadath_config_boost<BCO_NS_INFO> bconfig(std::string{fn});

  // get const EOS information - used for initializing EOS later
  const double h_cut = bconfig.eos<double>(HCUT);
  const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

  /* file containing KADATH fields must have same name as config file
   * with only the extension being different */
  std::string kadath_filename = bconfig.space_filename();

	FILE* fin = fopen(kadath_filename.c_str(), "r") ;
	Space_spheric_adapted space(fin) ;
	Scalar conf  (space, fin) ;
	Scalar lapse (space, fin) ;
  Vector shift (space, fin) ;
	Scalar logh  (space, fin) ;

	fclose(fin) ;

  // initialize container of references to scalar fields
  // for easy looping later
  std::vector<std::reference_wrapper<const Scalar>> quants;
  for (int i = 0; i < NUM_QUANTS; ++i)
    quants.push_back(std::cref(conf));

  quants[PSI] = std::cref(conf);
  quants[ALP] = std::cref(lapse);

  Base_tensor basis(shift.get_basis());
  Metric_flat fmet(space, basis);
  int ndom = space.get_nbr_domains();

  // fields depending on the coords
  CoordFields<Space_spheric_adapted> cf_generator(space);
  vec_ary_t coord_vectors {default_co_vector_ary(space)};

  // get origin of the system and initialize coordinate fields
  double xo = bco_utils::get_center(space,0);
  update_fields_co(cf_generator, coord_vectors, {}, xo);

  System_of_eqs syst(space, 0, ndom - 1);
  fmet.set_system(syst, "f");

  Param p;
  // add EOS user defined OPEs based on EOS type
  if(eos_type == "Cold_Table") {
    using namespace Kadath::Margherita;
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
                            2000 : bconfig.eos<int>(INTERP_PTS);

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
    syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
  }

  if(eos_type == "Cold_PWPoly") {
    using namespace Kadath::Margherita;
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
    syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
  } // end adding EOS OPEs

  // constants
  syst.add_cst("4piG", bconfig(BCO_QPIG));
  syst.add_cst("PI"  , M_PI);

  // coordinate dependent fields
  syst.add_cst("mg"  , *coord_vectors[GLOBAL_ROT]);
  syst.add_cst("ex"  , *coord_vectors[EX]);
  syst.add_cst("ey"  , *coord_vectors[EX]);
  syst.add_cst("ez"  , *coord_vectors[EX]);
  syst.add_cst("einf", *coord_vectors[S_INF]);
  syst.add_cst("sm",   *coord_vectors[S_BCO1]);

  // baryonic mass, dimensionless spin, central log enthalpy, angular frequency parameter
  syst.add_cst("Mb"  , bconfig(MB));
  syst.add_cst("chi" , bconfig(CHI));
  syst.add_cst("ome" , bconfig(OMEGA));
  syst.add_cst("Madm", bconfig(MADM));

  // gravitational fields
  syst.add_cst("P"   , conf);
  syst.add_cst("N"   , lapse);
  syst.add_cst("bet" , shift);

  // matter represented by log enthalpy
  syst.add_cst("H"   , logh);

  // common definitions
  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");

  // beta combined with the rotational part
  syst.add_def("omega^i = bet^i + ome * mg^i");

  // the extrinsic curvature
  syst.add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");
  
  syst.add_def("h = exp(H)");

  // definitions for the fluid 3-velocity
  syst.add_def("U^i = omega^i / N");

  Tensor A(syst.give_val_def("A"));
  Index ind(A);

  quants[AXX] = std::cref(A(ind));
  ind.inc();
  quants[AXY] = std::cref(A(ind));
  ind.inc();
  quants[AXZ] = std::cref(A(ind));
  ind.inc();
  ind.inc();
  quants[AYY] = std::cref(A(ind));
  ind.inc();
  quants[AYZ] = std::cref(A(ind));
  ind.inc();
  ind.inc();
  ind.inc();
  quants[AZZ] = std::cref(A(ind));

  quants[H] = std::cref(logh);

  Vector vel_kad(syst.give_val_def("U"));

  quants[UX] = std::cref(vel_kad(1));
  quants[UY] = std::cref(vel_kad(2));
  quants[UZ] = std::cref(vel_kad(3));

  quants[BETX] = std::cref(shift(1));
  quants[BETY] = std::cref(shift(2));
  quants[BETZ] = std::cref(shift(3));

  std::array<std::vector<double>,NUM_OUT> out;
  for(auto& v : out)
    v.resize(npoints);

  for (int i = 0; i < npoints; ++i) {
    std::vector<double> quant_vals(NUM_QUANTS);

    // get number of dimensions and construct index
    int ndim = 3;

    // construct Kadath point, shifted wrt the COM
    Point abs_coords(ndim);
    abs_coords.set(1) = xx[i];
    abs_coords.set(2) = yy[i];
    abs_coords.set(3) = zz[i];

    for (int k = 0; k < NUM_QUANTS; ++k) {
      quant_vals[k] = quants[k].get().val_point(abs_coords);
    }

    auto const psi = quant_vals[PSI];
    auto const psi2 = psi * psi;
    auto const psi4 = psi2 * psi2;

    out[ALPHA][i] = quant_vals[ALP];

    out[BETAX][i] = quant_vals[BETX];
    out[BETAY][i] = quant_vals[BETY];
    out[BETAZ][i] = quant_vals[BETZ];

    double g[3][3];
    g[0][0] = psi4;
    g[0][1] = 0.0;
    g[0][2] = 0.0;
    g[1][1] = psi4;
    g[1][2] = 0.0;
    g[2][2] = psi4;
    g[1][0] = g[0][1];
    g[2][0] = g[0][2];
    g[2][1] = g[1][2];

    out[GXX][i] = g[0][0];
    out[GXY][i] = g[0][1];
    out[GXZ][i] = g[0][2];
    out[GYY][i] = g[1][1];
    out[GYZ][i] = g[1][2];
    out[GZZ][i] = g[2][2];

    out[KXX][i] = quant_vals[AXX] * psi4;
    out[KXY][i] = quant_vals[AXY] * psi4;
    out[KXZ][i] = quant_vals[AXZ] * psi4;
    out[KYY][i] = quant_vals[AYY] * psi4;
    out[KYZ][i] = quant_vals[AYZ] * psi4;
    out[KZZ][i] = quant_vals[AZZ] * psi4;

    double h = std::exp(quant_vals[H]);
    double rho, eps, press;

    // get quantities point-wise, since h is smoothest, and cut data at h=1
    if(h == 1.) {
      rho = 0.;
      eps = 0.;
      press = 0.;
    }
    else {
      if(eos_type == "Cold_Table") {
        using namespace Kadath::Margherita;

        out[RHO][i] = EOS<Cold_Table, DENSITY>::get(h);
        out[EPS][i] = EOS<Cold_Table, EPSILON>::get(h);
        out[PRESS][i] = EOS<Cold_Table, PRESSURE>::get(h);
      }

      if(eos_type == "Cold_PWPoly") {
        using namespace Kadath::Margherita;

        out[RHO][i] = EOS<Cold_PWPoly, DENSITY>::get(h);
        out[EPS][i] = EOS<Cold_PWPoly, EPSILON>::get(h);
        out[PRESS][i] = EOS<Cold_PWPoly, PRESSURE>::get(h);
      }
    }

    out[VELX][i] = quant_vals[UX];
    out[VELY][i] = quant_vals[UY];
    out[VELZ][i] = quant_vals[UZ];
  } // for i
  return out;
}
