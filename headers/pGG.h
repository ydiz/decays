#pragma once

#include <qlat/qlat.h>
#include "io.h"
// #include "constants_macro.h"
#include "utils.h"


namespace Grid {
namespace QCD {


// to use Scidac, we cannot make PGGElem have only one level (like iMatrix<Complex, 4>)
// using PGGElem = iScalar<iScalar<iMatrix<Complex, 4>>>;
// using vPGGElem = iScalar<iScalar<iMatrix<vComplex, 4>>>;
// using LatticePGG = Lattice<vPGGElem>;

using PGGElem = iMatrix<iScalar<iScalar<Complex>>, 4>; // FIXME: I changed LatticePGG to be the same as LatticeKGG. Can I still use the program for pion????
using vPGGElem = iMatrix<iScalar<iScalar<vComplex>>, 4>;
using LatticePGG = Lattice<vPGGElem>;

using LatticeKGG = Lattice<iMatrix<iScalar<iScalar<vComplex> >, 4>>;
using LatticeLorentzColour = Lattice<iMatrix<iScalar<iMatrix<vComplex, 3> >, 4>>;

using LatticeFermionSite = typename LatticeFermion::vector_object::scalar_object;
using LatticePropagatorSite = typename LatticePropagator::vector_object::scalar_object;
using LatticeColourMatrixSite = typename LatticeColourMatrix::vector_object::scalar_object;
using LatticeComplexSite = typename LatticeComplex::vector_object::scalar_object;
using LatticePGGSite = typename LatticePGG::vector_object::scalar_object;
using LatticeKGGSite = typename LatticeKGG::vector_object::scalar_object;




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
