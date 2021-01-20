// Note: for free field test, the boundary condition must be periodic, i.e. boundary_phases must be (1,1,1,1), instead of (1,1,1,-1)

#include "../../../kaon/kaon.h"

// For 24ID, on 16 nodes, without using eigenvectors, need ~8h for 1 trajectory

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

#ifndef USE_CPS 
template<class FieldD,class FieldF, typename std::enable_if< getPrecision<FieldD>::value == 2, int>::type = 0,typename std::enable_if< getPrecision    <FieldF>::value == 1, int>::type = 0 > 
class MixedPrecisionConjugateGradientOp : public MixedPrecisionConjugateGradient<FieldD,FieldF>, public OperatorFunction<FieldD> { 
  public: 
    MixedPrecisionConjugateGradientOp(RealD tol, Integer maxinnerit, Integer maxouterit,  
        GridBase* _sp_grid, LinearOperatorBase<FieldF> &_Linop_f,  
        LinearOperatorBase<FieldD> &_Linop_d) : 
      MixedPrecisionConjugateGradient<FieldD,FieldF> (tol, maxinnerit, maxouterit, _sp_grid, _Linop_f, _Linop_d) {}; 

    void operator() (LinearOperatorBase<FieldD> &Linop, const FieldD &in, FieldD &out){ 
      this->MixedPrecisionConjugateGradient<FieldD,FieldF>::operator()(in,out); 
    } 
}; 
#endif



int main(int argc, char **argv)
{
  Grid_init(&argc,&argv);
  // zyd_init_Grid_Qlattice(argc, argv);

  // int Ls = 24;
  // double mass = 0.085, b = 2.5, M5 = 1.8;  // !! mass must be the mass of strange quark
  // vector<int> fdims = {24, 24, 24, 64};

  int Ls = 16;
  double mass = 0.04, b = 1.0, M5 = 1.0;   // light quark
  // double mass = 0.06, b = 1.0, M5 = 1.0;  // heavy quark
  vector<int> fdims = {8, 8, 8, 8};

  std::cout << "Lat size: " << fdims << std::endl;
  std::cout << "Ls: " << Ls << std::endl;
  std::cout << "M5: " << M5 << std::endl;
  std::cout << "mass: " << mass << std::endl;
  std::cout << "b: " << b << std::endl;

  string output_prefix = "./free_field_props";


  ////////////////////////////////////////////////////////////

  int T = fdims[3];

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(fdims, GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian *UGrid_f = SpaceTimeGrid::makeFourDimGrid(fdims, GridDefaultSimd(4, vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);

  GridCartesian *FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);
  GridCartesian *FGrid_f = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid_f);
  GridRedBlackCartesian *FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid_f);

  /////////////////////////////////////////
  // Read Gauge Field and do Gauge fixing
  /////////////////////////////////////////
  LatticeGaugeField Umu(UGrid);
  Umu = 1.0;   // Free field test
  // FieldMetaData header;
  // NerscIO::readConfiguration(Umu, header, Umu_dir + "/ckpoint_lat." + to_string(traj));


  /////////////////////////////////////////
  // Set Operators
  /////////////////////////////////////////

  LatticeGaugeFieldF Umu_f(UGrid_f);
  precisionChange(Umu_f, Umu);

  typename MobiusFermionD::ImplParams params;
  std::vector<Complex> boundary_phases(4, 1.);  // For free field test, all directions must be periodic
  // boundary_phases[3] = -1.;   // 
  params.boundary_phases = boundary_phases;

  MobiusFermionD Dmob = MobiusFermionD(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, b, b-1., params);
  MobiusFermionF Dmob_f = MobiusFermionF(Umu_f, *FGrid_f, *FrbGrid_f, *UGrid_f, *UrbGrid_f, mass, M5, b, b-1., params);
  SchurDiagTwoOperator<MobiusFermionD, LatticeFermionD> mobHermOp = SchurDiagTwoOperator<MobiusFermionD, LatticeFermionD>(Dmob);
  SchurDiagTwoOperator<MobiusFermionF, LatticeFermionF> mobHermOp_f = SchurDiagTwoOperator<MobiusFermionF, LatticeFermionF>(Dmob_f);

  double resid = 1e-8;
  int max_iters = 10000;
  MixedPrecisionConjugateGradientOp<LatticeFermionD, LatticeFermionF> mCG(resid, max_iters, 50, FrbGrid_f, mobHermOp_f, mobHermOp);
  SchurRedBlackDiagTwoSolve<LatticeFermionD> solver(mCG);

  /////////////////////////////////////////
  // Calculate Wall Source Propagators
  /////////////////////////////////////////

  Lattice<iScalar<vInteger>> t(UGrid);
  LatticeCoordinate(t, Tp);

  int tW = 0; // position of source is tW=0
  // for(int tW=0; tW<T; ++tW) {
  std::cout << GridLogMessage << "t_wall: " << tW << std::endl;
  LatticePropagator fullSrc(UGrid);
  fullSrc = 1.;
  fullSrc = where((t == tW), fullSrc, 0. * fullSrc);

  LatticePropagator prop(UGrid);
  for(int s = 0; s < 4; ++s) {
    for(int c = 0; c < 3; ++c) {
      std::cout << GridLogMessage << "s: " << s << " c: " << c << std::endl;
      LatticeFermion src(UGrid), src5d(FGrid), sol(UGrid), sol5d(FGrid);

      PropToFerm<WilsonImplR>(src, fullSrc, s, c);

      Dmob.ImportPhysicalFermionSource(src, src5d);

      solver(Dmob, src5d, sol5d); // zyd: Dmob is not used, but required for syntax // In the MixedPrecisionConjugateGradientOp class, Dmob is not used

      Dmob.ExportPhysicalFermionSolution(sol5d, sol); // 5d->4d, v(x) = P_L v(x,0) + P_R v(x, Ls-1)

      FermToProp<WilsonImplR>(prop, sol, s, c);
    }
  }

  // string prop_fname = output_prefix + "/wall_s/" + std::to_string(tW);
  string prop_fname = output_prefix + "/wall_l/" + std::to_string(tW);
  writeScidac_prop_d2f(prop, prop_fname);

  // // std::cout << prop << std::endl;
  // std::cout << "Full propagator at {0,0,0,0}" << std::endl;
  // print_grid_field_site(prop, {0,0,0,0});


  LatticeSpinMatrix tmp = peekColour(prop, 0, 0); // The color index dependence of prop is delta_{a,b}, and is trivial
  for(int t=0; t<T; ++t) print_grid_field_site(tmp, {0,0,0,t}); // I checked that all sites with the same t are equal

  // }  // end of loop of tW


  std::cout << "FINISHED!" << std::endl;

  Grid_finalize();
}

