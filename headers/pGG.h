#pragma once

#include <qlat/qlat.h>
#include "io.h"
// #include "constants_macro.h"
#include "utils.h"


namespace Grid {
namespace QCD {

// to use Scidac, we cannot make PGGElem have only one level (like iMatrix<Complex, 4>)
using PGGElem = iScalar<iScalar<iMatrix<Complex, 4>>>;
using vPGGElem = iScalar<iScalar<iMatrix<vComplex, 4>>>;
using LatticePGG = Lattice<vPGGElem>;

// translational factor
void get_translational_factor(LatticeComplex &lat, double Mpi) {

	parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++){
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

		double val;
		int xt = qlat::smod(gcoor[Tdir], lat.Grid()->_fdimensions[Tdir]);
    val = std::exp( 0.5 * Mpi * xt); // translation factor

		typename LatticeComplex::vector_object::scalar_object m;
		m()()() = Complex(val, 0.);
		pokeLocalSite(m, lat, lcoor);
	}


}


}}
