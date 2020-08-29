#include "../kaon/kaon.h"

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
  // Grid_init(&argc,&argv);
  zyd_init_Grid_Qlattice(argc, argv);

  int Ls = 24;
  double mass = 0.085, b = 2.5, M5 = 1.8;  // !! mass must be the mass of strange quark
  vector<int> fdims = {24, 24, 24, 64};
  // string output_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/propagators/24ID";
  string output_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_my_props";
  string Umu_dir = "/hpcgpfs01/work/lqcd/etap/chulwoo/evec/24ID1Gev/configurations";
  string gt_dir = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/gauge_transform";  // gauge transformations


  int target_traj;
  if( GridCmdOptionExists(argv, argv+argc, "--traj") ) {
    string arg = GridCmdOptionPayload(argv, argv+argc, "--traj");
    GridCmdOptionInt(arg, target_traj);
  }
  else {
    std::cout << "traj not specified; exiting" << std::endl;
    assert(0);
  }

  int traj_start = target_traj;
  int traj_end = target_traj;
  int traj_sep = 10;

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

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    /////////////////////////////////////////
    // Read Gauge Field and do Gauge fixing
    /////////////////////////////////////////
    LatticeGaugeField Umu(UGrid);
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, Umu_dir + "/ckpoint_lat." + to_string(traj));

    ////////////////////////////////////////
    // Option 1. Read Luchang's gauge transformation
    ////////////////////////////////////////
    // string gt_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/propagators/results=" + to_string(traj) + "/gauge-transform.field";
    string gt_path = gt_dir + "/gauge-transform.field." + to_string(traj);
    LatticeColourMatrix gt(UGrid);
    read_luchang_gt(gt, gt_path);
    SU<3>::GaugeTransform(Umu, gt);  // I checked, gt(x) U gt(x+mu)^dagger  is in Coulomb gauge; (By running GaugeFixer, the Phi an Omega are small)
    
    ////////////////////////////////////////
    // Option 2. Calculate gauge transformation
    // Do not converge!!!! No matter whether I turn off Fourier or not.  
    // Even if I start with the configuration that has already been gauge fixed, it will stray away.
    // Maybe I need to use smaller alpha (alpha is step size)
    ////////////////////////////////////////
    // double alpha = 0.1;
    // int coulomb_dir = Nd - 1;
    // LatticeColourMatrix tmp_gt(UGrid);
    // FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu, tmp_gt, alpha, 10000, 1.0e-12, 1.0e-12, false, coulomb_dir); // Note: it is very possible that gauge fixing does not converge; to converge, Phi and Omega must be smaller than their tolerance.
    //
    // string gt_fname = output_prefix + "/gt/" + to_string(traj);
    // writeScidac(gt, gt_fname);

    ///////////////////////////////////////////////


    /////////////////////////////////////////
    // Set Operators
    /////////////////////////////////////////

    LatticeGaugeFieldF Umu_f(UGrid_f);
    precisionChange(Umu_f, Umu);

    // Must set boundary phase in time direction to -1
    typename MobiusFermionD::ImplParams params;
    std::vector<Complex> boundary_phases(4, 1.);
    boundary_phases[3] = -1.;
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

    for(int tW=0; tW<T; ++tW) {
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

      string prop_fname = output_prefix + "/wall_s/" + to_string(traj)  + "/" + std::to_string(tW);
      writeScidac_prop_d2f(prop, prop_fname);
    }  // end of loop of tW

  }  // end of loop of traj


  Grid_finalize();
}

