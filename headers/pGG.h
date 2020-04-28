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

using LatticePropagatorSite = typename LatticePropagator::vector_object::scalar_object;
using LatticeComplexSite = typename LatticeComplex::vector_object::scalar_object;
using LatticePGGSite = typename LatticePGG::vector_object::scalar_object;




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

// translational factor
void get_luchang_exp_factor(LatticeComplex &lat, double Mpi, double t_min) {

	parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++){
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

		double val;
		int xt = qlat::smod(gcoor[Tdir], lat.Grid()->_fdimensions[Tdir]);

    if(xt>=0) val = std::exp( Mpi * t_min);
    else val = std::exp( Mpi * (t_min-xt) );
    // val = std::exp( 0.5 * Mpi * xt); // translation factor

		typename LatticeComplex::vector_object::scalar_object m;
		m()()() = Complex(val, 0.);
		pokeLocalSite(m, lat, lcoor);
	}


}





}}
