#pragma once

#include <headers/headers.h>

namespace Grid {
namespace QCD {

void form_factor_integrand(LatticePGG &lat, double Mpi_lat) {

	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){

    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    gcoor = my_smod(gcoor, lat._grid->_fdimensions);

    double w = std::sqrt(gcoor[0]*gcoor[0] + gcoor[1]*gcoor[1] + gcoor[2]*gcoor[2]);

    double val;
    if(w==0) val = 0.;
    else {
      double t = Mpi_lat * 0.5 * w;
      double sin_t, cos_t;
      sincos(t, &sin_t, &cos_t);
      val = ( - cos_t * Mpi_lat * w + 2 * sin_t) / w / w / w * std::exp(Mpi_lat * 0.5 * gcoor[Tdir]);
    }

		typename LatticePGG::vector_object::scalar_object m;
		m = 0.;
		m()()(0, 1) = Complex(val * gcoor[Zdir], 0); 
		m()()(0, 2) = Complex(-val * gcoor[Ydir], 0); // Minus sign comes from spinor matrix
		m()()(1, 2) = Complex(val * gcoor[Xdir], 0); 
		m()()(1, 0) = - m()()(0, 1);
		m()()(2, 0) = - m()()(0, 2);
		m()()(2, 1) = - m()()(1, 2);

		pokeLocalSite(m, lat, lcoor);
	}
  // std::cout << lat << std::endl;
  // exit(0);
}



std::vector<double> form_factor(const LatticePGG &three_point, const LatticePGG &leptonic, double hadron_coeff, double lep_coeff) {

	std::vector<double> ret = mult_HL_cutoff(three_point, leptonic);

  // double other_coeff = 8 * M_PI * M_PI * Fpi / Mpi_lat / Mpi_lat / Mpi_lat / Mpi_lat;
  std::vector<double> F(ret.size());
  for(int i=0; i<ret.size(); ++i) F[i] = hadron_coeff * lep_coeff * ret[i];

  return F;
}



}}

