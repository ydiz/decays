#pragma once

#include "utils.h"
#include "amplitude.h"


namespace Grid {
namespace QCD {

void form_factor_integrand_scalar(LatticeComplex &lat) {

	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){

    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    gcoor = my_smod(gcoor, lat._grid->_fdimensions);

    double w = std::sqrt(gcoor[0]*gcoor[0] + gcoor[1]*gcoor[1] + gcoor[2]*gcoor[2]);

    double val;
    if(w==0) val = 0.;
    else {
      double t = M_PION * 0.5 * w;
      double sin_t, cos_t;
      sincos(t, &sin_t, &cos_t);
      val = ( - cos_t * M_PION * w + 2 * sin_t) / w / w / w * std::exp(M_PION * 0.5 * gcoor[Tdir]);
    }

		typename LatticeComplex::vector_object::scalar_object m;
		m()()() = Complex(val, 0); 

		pokeLocalSite(m, lat, lcoor);
	}
}



LatticeComplex toScalar(const LatticePGG& pgg) {
  LatticeComplex ret(pgg._grid);
	parallel_for(int ss=0; ss<pgg._grid->lSites(); ss++){

    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(pgg._grid, ss, lcoor, gcoor);

    gcoor = my_smod(gcoor, pgg._grid->_fdimensions);

		typename LatticePGG::vector_object::scalar_object pgg_site;
    peekLocalSite(pgg_site, pgg, lcoor);
		typename LatticeComplex::vector_object::scalar_object ret_site;
    ret_site()()() += double(gcoor[2]) * (pgg_site()()(0, 1) - pgg_site()()(1, 0)) - double(gcoor[1]) * (pgg_site()()(0,2) - pgg_site()()(2,0)) + double(gcoor[0]) * (pgg_site()()(1,2) - pgg_site()()(2,1));
		pokeLocalSite(ret_site, ret, lcoor);
  }

  return ret;

}

// std::vector<double> form_factor_scalar_H(const LatticePGG &three_point, const LatticePGG &leptonic) {
//   // static LatticeComplex pp(three_point._grid); 
//   // static bool pp_initialzed = false;
//   // if(!pp_initialzed) {
//   //   get_translational_factor(pp);// translational factor 
//   //   pp_initialzed = true;
//   // }
//
// 	LatticePGG hadronic(three_point._grid);
// 	// hadronic = three_point * pp; // Do not need to shift hadronic part
//   double Z_V = 0.72672;
//   double hadron_coeff = 1./ (3 * std::sqrt(2)) * Z_V * Z_V * 2 * M_PION / (9428.492 / std::sqrt(32*32*32.));
//   hadronic = imag(three_point); 
//   
//   LatticeComplex H(three_point._grid);
//   H = toScalar(hadronic);
//   H = hadron_coeff * H;
//   print_grid_field_site(H, std::vector<int>{1,8,13,5});
//
// 	std::vector<double> ret = mult_HL_cutoff(hadronic, leptonic);
//
// 	double me = 511000;
//   double other_coeff = 8 * M_PI * M_PI * 0.092424 / M_PION / M_PION / M_PION / M_PION;
//   std::vector<double> F(ret.size());
//   for(int i=0; i<ret.size(); ++i) F[i] = hadron_coeff * other_coeff * ret[i];
//
//   return F;
// }



}}

