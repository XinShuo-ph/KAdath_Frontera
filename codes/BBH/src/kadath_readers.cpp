#include "kadath_bin_bh.hpp"
#include "python_reader.hpp"
#include "coord_fields.hpp"
#include "Configurator/config_binary.hpp"
#include "bco_utilities.hpp"
#include <string>
typedef Kadath::Space_bin_bh space_t;
struct bbh_vars_t : public Kadath::vars_base_t<bbh_vars_t> {};

template<> Kadath::var_vector Kadath::vars_base_t<bbh_vars_t>::vars = {
                                                                        {"conf", SCALAR},
                                                                        {"lapse", SCALAR},
                                                                        {"shift", VECTOR},
                                                                      };


class bbh_reader_t : public Kadath::python_reader_t<space_t, bbh_vars_t> {
  std::string config_filename;
  kadath_config_boost<BIN_INFO> bconfig;
  public:
  bbh_reader_t(std::string const filename) : Kadath::python_reader_t<space_t, bbh_vars_t>(filename),
                                             config_filename(filename.substr(0,filename.size()-3)+"info"),
                                             bconfig(config_filename) {
    this->compute_defs();
  }

  void compute_defs() {
    using namespace Kadath;

    Kadath::Scalar const & conf  = extractField<Kadath::Scalar>("conf");
    Kadath::Scalar const & lapse = extractField<Kadath::Scalar>("lapse");
    Kadath::Vector const & shift = extractField<Kadath::Vector>("shift");

    Base_tensor basis(shift.get_basis());
  	Metric_flat fmet (space, basis) ;

  	int ndom = space.get_nbr_domains() ;

   	double xc1 = bco_utils::get_center(space,space.BH1);
   	double xc2 = bco_utils::get_center(space,space.BH2);
    double xo  = bco_utils::get_center(space,ndom-1);

  	// setup coordinate vector fields for System_of_eqs
    CoordFields<Space_bin_bh> cf_generator(space);
    std::array<Vector*, NUM_VECTORS> coord_vectors {};
    for(auto& el : coord_vectors) el = new Vector(space,CON,basis);
    update_fields(cf_generator, coord_vectors, {}, xo, xc1, xc2);
    
    Vector CART(space, CON, basis);
    CART = cf_generator.cart();
    // end coordinate field setup
    
  	// Setup system of equations and definitions
    System_of_eqs syst (space, 0, ndom-1) ;

  	fmet.set_system(syst, "f") ;
    syst.add_cst    ("PI"   , M_PI);

  	// coordinate fields
    syst.add_cst    ("mg"   , *coord_vectors[GLOBAL_ROT]) ;
  	syst.add_cst    ("ex"   , *coord_vectors[EX])  ;
  	syst.add_cst    ("ey"   , *coord_vectors[EY])  ;
  	syst.add_cst    ("ez"   , *coord_vectors[EZ])  ;
  	syst.add_cst    ("einf" , *coord_vectors[S_INF]) ;
  
  	// BCO centered fields
    syst.add_cst    ("mm"   , *coord_vectors[BCO1_ROT]) ;
  	syst.add_cst    ("mp"   , *coord_vectors[BCO2_ROT]) ;
  	syst.add_cst    ("sm"   , *coord_vectors[S_BCO1]) ;
  	syst.add_cst    ("sp"   , *coord_vectors[S_BCO2]) ;
  
	  // BH1 Quantities
    syst.add_cst    ("Mm"   , bconfig(MIRR,BCO1)) ;
	  syst.add_cst    ("chim" , bconfig(CHI,BCO1)) ;
    syst.add_cst    ("CMm"  , bconfig(MCH,BCO1));

	  // BH2 Quantities
    syst.add_cst    ("Mp"   , bconfig(MIRR,BCO2)) ;
	  syst.add_cst    ("chip" , bconfig(CHI,BCO2)) ;
    syst.add_cst    ("CMp"  , bconfig(MCH,BCO2));

    // Binary Quantities
    syst.add_cst    ("ome"  , bconfig(GOMEGA));
    syst.add_cst    ("xaxis", bconfig(COM));
    syst.add_cst    ("yaxis", bconfig(COMY));
    
    // Boosts are not commonly used so we check
    // if they need to be added
    std::string boosts{};
    if(!std::isnan(bconfig.set(BVELX))) {
      syst.add_cst    ("xvel" , bconfig(BVELX));
      boosts+=" + xvel * ex^i";
    }
    if(!std::isnan(bconfig.set(BVELY))) {
      syst.add_cst    ("yvel" , bconfig(BVELY));
      boosts+=" + yvel * ey^i";
    }

    // check for ADOT so we don't get errors.
    if(!std::isnan(bconfig.set(ADOT))) {
      syst.add_cst    ("adot" , bconfig(ADOT));
      syst.add_cst    ("r"   , CART);
      syst.add_def    ("comr^i = r^i - xaxis * ex^i + yaxis * ey^i");
      boosts+=" + adot * comr^i";
    }
   
    // Fields
    syst.add_cst    ("P"    , conf);
    syst.add_cst    ("N"    , lapse);
    syst.add_cst    ("bet"  , shift);
  
    // Helpful definitions
    syst.add_def    ("NP      = P*N");
    syst.add_def    ("Ntilde  = N / P^6");

    // Define COM shifted orbital rotation field
    syst.add_def    ("Morb^i  = mg^i + xaxis * ey^i + yaxis * ex^i");
    // Define global shift (inertial + corotating)
    std::string B { "B^i     = bet^i + ome * Morb^i" + boosts };
    syst.add_def    (B.c_str());

    // Conformal Extrinsic curvature.
    syst.add_def    ("A^ij    = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

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

    syst.add_def    ("Bx      = B^i  * ex_i ");
    syst.add_def    ("By      = B^i  * ey_i");
    syst.add_def    ("Bz      = B^i  * ez_i");
  
    // Local spin definition
    syst.add_def    ("intSm   = A_ij * mm^i * sm^j   / 8 / PI") ;
    syst.add_def    ("intSp   = A_ij * mp^i * sp^j   / 8 / PI") ;
  
    // Irreducible mass integrand
    syst.add_def    ("intMsq  = P^4 / 16. / PI") ;

    // ADM Quantities
    syst.add_def    (ndom - 1, "intJ = multr(A_ij * Morb^j * einf^i) / 8 / PI");
    syst.add_def    (ndom - 1, "Madm = -dr(P) / 2 / PI");
    syst.add_def    ("intPx   = A_ij * ex^j * einf^i / 8 / PI") ;
    syst.add_def    ("intPy   = A_ij * ey^j * einf^i / 8 / PI") ;
    syst.add_def    ("intPz   = A_ij * ez^j * einf^i / 8 / PI") ;

    // Komar mass measured at inf
    syst.add_def    (ndom - 1, "Mk   =  dr(N) / 4 / PI");
 
    // Constraint equations
    syst.add_def    ("eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8") ;
    syst.add_def    ("eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def    ("eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j Ntilde");
    // end system of equations

    // Populate Boost Dictionary
    double S1     = space.get_domain(space.BH1+2)->integ(syst.give_val_def("intSm")()(space.BH1+2) , OUTER_BC); 
    double S2     = space.get_domain(space.BH2+2)->integ(syst.give_val_def("intSp")()(space.BH2+2) , OUTER_BC);
    double Mirr1  = sqrt(space.get_domain(space.BH1+2)->integ(syst.give_val_def("intMsq")()(space.BH1+2) , INNER_BC));
    double Mirr2  = sqrt(space.get_domain(space.BH2+2)->integ(syst.give_val_def("intMsq")()(space.BH2+2) , INNER_BC));

    // Scalar fields of the constraint equations
    vars["cP"]    = syst.give_val_def("eqP");
    vars["cNP"]    = syst.give_val_def("eqNP");

    // computed quantities on the BHs
    vars["Mirr1"] = Mirr1;
    vars["Mirr2"] = Mirr2;
    vars["S1"]    = S1; 
    vars["S2"]    = S2;
    vars["Mch1"]  = std::sqrt(Mirr1 * Mirr1 + S1 * S1 / 4. / Mirr1 / Mirr1);
    vars["Mch2"]  = std::sqrt(Mirr2 * Mirr2 + S2 * S2 / 4. / Mirr2 / Mirr2);

    // ADM and Komar quantities
    vars["Jadm"]  = space.get_domain(ndom-1)->integ(syst.give_val_def("intJ")()(ndom-1) , OUTER_BC);
    vars["Madm"]  = space.get_domain(ndom-1)->integ(syst.give_val_def("Madm")()(ndom-1) , OUTER_BC);
    vars["Mk"]    = space.get_domain(ndom-1)->integ(syst.give_val_def("Mk")()(ndom-1) , OUTER_BC);
    vars["Px"]    = space.get_domain(ndom - 1)->integ(syst.give_val_def("intPx")()(ndom-1), OUTER_BC);
    vars["Py"]    = space.get_domain(ndom - 1)->integ(syst.give_val_def("intPy")()(ndom-1), OUTER_BC);
    vars["Pz"]    = space.get_domain(ndom - 1)->integ(syst.give_val_def("intPz")()(ndom-1), OUTER_BC);

    // BCO coordinate centers
    vars["xc1"]   = xc1;
    vars["xc2"]   = xc2;

    // Conformal Aij contractions
    vars["Axx"]   = syst.give_val_def("Axx");
    vars["Ayy"]   = syst.give_val_def("Ayy");
    vars["Azz"]   = syst.give_val_def("Azz");
    vars["Axy"]   = syst.give_val_def("Axy");
    vars["Ayz"]   = syst.give_val_def("Ayz");
    vars["Axz"]   = syst.give_val_def("Axz");
    vars["A"]     = syst.give_val_def("TraceA");

    // inertial shift contractions
    vars["betx"]  = syst.give_val_def("betx");
    vars["bety"]  = syst.give_val_def("bety");
    vars["betz"]  = syst.give_val_def("betz");

    // total shift contractions
    vars["Bx"]    = syst.give_val_def("Bx");
    vars["By"]    = syst.give_val_def("By");
    vars["Bz"]    = syst.give_val_def("Bz");

    //create old radius scalar field
    Scalar space_radius(space);
    space_radius.annule_hard();
    for(int d = 0; d < ndom; ++d){
      space_radius.set_domain(d) = space.get_domain(d)->get_radius();
    }
    space_radius.std_base();

    auto find_rmin = [&](int d) {
			Index pos(space.get_domain(d)->get_nbr_points());
		  pos.set_start();

  		double rmin = space_radius(d)(pos);
  		do {
				pos.set(0) = space.get_domain(d)->get_nbr_points()(0) - 1;
				double r = space_radius(d)(pos);

				if(r < rmin)
					rmin = r;
		  } while(pos.inc());
		  return rmin;

    };

    // create a list of domain radii
    boost::python::list d_outer_r;
    for(int d = 0; d < ndom; ++d)
			d_outer_r.append(find_rmin(d));

		vars["dom_r"] = d_outer_r;

    // This exports the radii only related to specific domains.  In this case
    // we care only about the domains in the vicinity of a BH.
    // For example, this is helpful when looking at constraint violations
    // near a BH and determining if the domains are having a positive or negative
    // impact on the solution especially when additional shells are added
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
    export_radii(space.BH1, space.BH2, "BH1_R");
    export_radii(space.BH2, space.OUTER, "BH2_R");
    
    // this was used to do a colored plot of the domains
    // like kids coloring books - color the numbered areas.
    Scalar dom_colors(conf);
    dom_colors.annule_hard();
    int c = 1;
    for(int d = 0; d < ndom; ++d){
        //set BH interiors to the same color
        if( d == space.BH1    ||
            d == space.BH1+1  ||
            d == space.BH2    ||
            d == space.BH2+1){
          dom_colors.set_domain(d) = 0;
        }
        //set inner_adapted domains to the same color
        else if( d == space.BH1+2 || d == space.BH2+2){
          dom_colors.set_domain(d) = 1;
        }
        //set chi_first domains to the same color
        else if( d == space.OUTER || d == space.OUTER+4) {
          dom_colors.set_domain(d) = 2;
        }
        //rect domains
        else if( d == space.OUTER+1 || d == space.OUTER+3){
          dom_colors.set_domain(d) = 3;
        }
        //eta + shells + compactified domains
        else {
          dom_colors.set_domain(d) = 3 + c;
          c++;
        }
    }
    dom_colors.std_base();
    vars["dom_color_chart"] = dom_colors;
    vars["rfield"] = space_radius;
  // clean-up vector fields
  for(auto& el : coord_vectors) delete el;
 }
};

BOOST_PYTHON_MODULE(kadath_readers)
{
    // initialize python types
    Kadath::initPythonBinding<space_t>();
    Kadath::constructPythonReader<bbh_reader_t>("bbh_reader");
}
