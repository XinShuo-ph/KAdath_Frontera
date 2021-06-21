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
#include "kadath_adapted_polar.hpp"
#include "kadath_adapted.hpp"
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"
#include "coord_fields.hpp"
#include "mpi.h"

using namespace Kadath;

template<class eos_t, typename config_t>
void reader_3d(config_t bconfig);

int main(int argc, char **argv) {
  // initialize MPI
  int rc = MPI_Init(&argc, &argv);

  if (rc != MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  // expecting a configuration file on execution
  if(argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<ID base name>.info" << std::endl;
    std::cerr << "e.g. ./reader converged.NS.9.info" << endl;
    std::_Exit(EXIT_FAILURE);
  }

  // load the configuration
  std::string ifilename{argv[1]};
  kadath_config_boost<BCO_NS_INFO> bconfig(ifilename);

  // setup the EOS
  const double h_cut = bconfig.eos<double>(HC);
  const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

  if(eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut);

    // call reader to output diagnostics
    if(bconfig(DIM) == 3) reader_3d<eos_t>(bconfig);
  } else if(eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
                            2000 : bconfig.eos<int>(INTERP_PTS);

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);

    // call reader to output diagnostics
    if(bconfig(DIM) == 3) reader_3d<eos_t>(bconfig);
  } else {
    std::cerr << "Unknown EOSTYPE." << endl;
    std::_Exit(EXIT_FAILURE);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}


template<class eos_t, typename config_t>
void reader_3d(config_t bconfig) {
  // reconstruction of the data filename
  std::string ifilename{bconfig.config_filename()};
  int idx = ifilename.rfind(".");
  std::string spacein = bconfig.config_outputdir()+'/'+ifilename.substr(0,idx)+".dat";

  // load domain decomposition and fields from binary file
  FILE* ff1 = fopen(spacein.c_str(), "r") ;
	Space_spheric_adapted space (ff1) ;
	Scalar conf   (space, ff1) ;
	Scalar lapse  (space, ff1) ;
  Vector shift  (space, ff1) ;
  Scalar logh   (space, ff1) ;
	fclose(ff1) ;

  // number of dimensions, a cartesian type of basis and the flat bg metric  
  int ndom = space.get_nbr_domains();
  Base_tensor basis(space, CARTESIAN_BASIS);
  Metric_flat fmet(space, basis);

  // fields depending on the coords
  CoordFields<Space_spheric_adapted> cf_generator(space);
  std::array<Vector*, NUM_VECTORS> coord_vectors {};
  coord_vectors[GLOBAL_ROT] = new Vector(space,CON,basis);
  coord_vectors[EX]         = new Vector(space,CON,basis);
  coord_vectors[EY]         = new Vector(space,CON,basis);
  coord_vectors[EZ]         = new Vector(space,CON,basis);
  coord_vectors[S_BCO1]     = new Vector(space,CON,basis);
  coord_vectors[S_INF]      = new Vector(space,CON,basis);

  // get origin of the system and initialize coordinate fields
  double xo = bco_utils::get_center(space,0);
  update_fields_co(cf_generator, coord_vectors, {}, xo);

  // central values of the matter fields
  double loghc = bco_utils::get_boundary_val(0, logh, INNER_BC);
  double hc = std::exp(loghc);
  double nc = EOS<eos_t,DENSITY>::get(hc);
  double pc = EOS<eos_t,PRESSURE>::get(hc);

  // minimal and maximal radius of the adapted surface domain
  auto [ rmin, rmax ] = bco_utils::get_rmin_rmax(space, 1);
  // inner radius of the nucleus
  double rin1 = bco_utils::get_radius(space.get_domain(0), OUTER_BC);

  // setup the system of equations to define derived quantities 
  System_of_eqs syst(space, 0, ndom - 1);

  // flat background metric
  fmet.set_system(syst, "f");

  // constants
  syst.add_cst("4piG", bconfig(BCO_QPIG));
  syst.add_cst("PI"  , M_PI);

  // baryonic mass, dimensionless spin, central log enthalpy, angular frequency parameter
  syst.add_cst("Mb"  , bconfig(MB));
  syst.add_cst("chi" , bconfig(CHI));
  syst.add_cst("Hc"  , loghc);
  syst.add_cst("ome" , bconfig(OMEGA));
  syst.add_cst("Madm", bconfig(MADM));

  // coordinate dependent fields
  syst.add_cst("mg"  , *coord_vectors[GLOBAL_ROT]);
  syst.add_cst("ex"  , *coord_vectors[EX]);
  syst.add_cst("ey"  , *coord_vectors[EX]);
  syst.add_cst("ez"  , *coord_vectors[EX]);
  syst.add_cst("einf", *coord_vectors[S_INF]);

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

  // extrinsic curvature
  syst.add_def("A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

  // integrants at infinity
  syst.add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 8. / PI");
  syst.add_def(ndom - 1, "intMadmalt = -dr(P) / 2 / PI");
  syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P / 2 / PI");
  syst.add_def(ndom - 1, "intMk = (einf^i * D_i N - A_ij * einf^i * bet^j) / 4 / PI");

  // setup operators for the EOS
  Param p;
  syst.add_ope ("eps"   , &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope ("pres" , &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope ("rho"   , &EOS<eos_t,DENSITY>::action, &p);

  // setup derived matter quantities
  for (int d = 0; d < ndom; d++) {
    switch (d) {
    case 0:
    case 1:
      syst.add_def(d, "h = exp(H)");

      syst.add_def(d, "U^i = omega^i / N");
      syst.add_def(d, "Usquare = P^4 * U_i * U^i");
      syst.add_def(d, "Wsquare = 1. / (1. - Usquare)");
      syst.add_def(d, "W = sqrt(Wsquare)");

      syst.add_def(d, "intMb = P^6 * rho(h) * W");
      syst.add_def(d, "firstint = H + log(N) - log(W)");
      break;
    }
  }

  // ADM angular momentum at infinity 
  Val_domain integJ(syst.give_val_def("intJ")()(ndom - 1));
  double J = space.get_domain(ndom - 1)->integ(integJ, OUTER_BC);

  // two (equivalent) computations of the ADM mass at infinity
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  Val_domain integMadmalt(syst.give_val_def("intMadmalt")()(ndom - 1));
  double Madmalt = space.get_domain(ndom - 1)->integ(integMadmalt, OUTER_BC);

  // Komar mass at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  // baryonic mass as volume integral over the two stellar domains 
  double baryonic_mass = syst.give_val_def("intMb")()(0).integ_volume() +
                         syst.give_val_def("intMb")()(1).integ_volume();

  // integrate proper area to get areal radius
  // this is only true for constant radius fields
  // for rotating stars, this is not accurate and was only used
  // as a diagnostic for irrotation stars.
  double AR = 0.;
  if(J <= 1e-12) {
    syst.add_def ("intMsq = P^4");
    double A = space.get_domain(2)->integ(syst.give_val_def("intMsq")()(2), INNER_BC);
    AR = sqrt(A / 4. / acos(-1.));
  }

  // output to stdout
  std::cout << bconfig;
  std::cout << "******************************************\n";
  #define FORMAT std::setw(20) << std::right << std::setprecision(5) << std::fixed << std::showpos
  std::cout << FORMAT << "RES = "             << space.get_domain(0)->get_nbr_points()(0)   << std::endl
            << FORMAT << "Nucleus Radius = "  << rin1                                       << std::endl
            << FORMAT << "Surface Radii = {"  << rmin << ", " << rmax << "}"                << std::endl;
  if(AR > 0.)
    std::cout << FORMAT << "Areal Radius = "    << AR                                         << std::endl;
  std::cout << FORMAT << "HC = "              << hc                                         << std::endl
            << FORMAT << std::setprecision (8)<< std::scientific << "NC = " << nc           << std::endl
            << FORMAT << std::setprecision (8)<< std::scientific << "PC = " << pc           << std::endl
            << FORMAT << "Madm = "            << Madm                                       << std::endl
            << FORMAT << "Mk = "              << Mk                                         << std::endl
            << FORMAT << "Mb = "              << baryonic_mass                              << std::endl
            << FORMAT << "Jadm = "            << J                                          << std::endl;
  
  // cleanup coord field pointers
  for(auto& el : coord_vectors) delete el;
}
