// p.s. To switch to heavy quark, need to change 1) mass parameter 2) prop_name (saved parameter path)
//
// Note: for free field test, the boundary condition must be periodic, i.e. boundary_phases must be (1,1,1,1), instead of (1,1,1,-1)

#include "../../../kaon/kaon.h"

using namespace std;
using namespace Grid;

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

  double mass = 0.04, b = 1.0, M5 = 1.0; // light quark
  string quark_type = "l";  // used in output path
  // double mass = 0.06, b = 1.0, M5 = 1.0; // heavy quark
  // string quark_type = "s";  // used in output path
  int Ls = 16;
  vector<int> fdims = {8, 8, 8, 8};
  // double time_boundary_phase = 0.; // No boundary phase; periodic boundary condition
  double time_boundary_phase = 113. * M_PI / 180.;  // 113 degrees

  std::cout << "Lat size: " << fdims << std::endl;
  std::cout << "Ls: " << Ls << std::endl;
  std::cout << "M5: " << M5 << std::endl;
  std::cout << "mass: " << mass << std::endl;
  std::cout << "quark_type: " << quark_type << std::endl;
  std::cout << "b: " << b << std::endl;
  std::cout << "time_boundary_phase: " << time_boundary_phase << std::endl;

  // string output_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_my_props";
  string output_prefix = "./free_field_props";

  ////////////////////////////////////////////////////////////

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

  // Must set boundary phase in time direction to -1
  typename MobiusFermionD::ImplParams params;
  std::vector<Complex> boundary_phases(4, 1.);   // For free field test, all directions must be periodic
  boundary_phases[3] = std::exp(Complex(0., time_boundary_phase)); // exp(i alpha)
  std::cout << "boundary_phases: " << boundary_phases << std::endl;
  // boundary_phases[3] = -1.;
  params.boundary_phases = boundary_phases;

  MobiusFermionD Dmob = MobiusFermionD(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, b, b-1., params);
  MobiusFermionF Dmob_f = MobiusFermionF(Umu_f, *FGrid_f, *FrbGrid_f, *UGrid_f, *UrbGrid_f, mass, M5, b, b-1., params);
  SchurDiagTwoOperator<MobiusFermionD, LatticeFermionD> mobHermOp = SchurDiagTwoOperator<MobiusFermionD, LatticeFermionD>(Dmob);
  SchurDiagTwoOperator<MobiusFermionF, LatticeFermionF> mobHermOp_f = SchurDiagTwoOperator<MobiusFermionF, LatticeFermionF>(Dmob_f);

  // double resid = 1e-8;   // FIXME: change to 1e-16
  double resid = 1e-16;
  int max_iters = 10000;
  MixedPrecisionConjugateGradientOp<LatticeFermionD, LatticeFermionF> mCG(resid, max_iters, 50, FrbGrid_f, mobHermOp_f, mobHermOp);
  SchurRedBlackDiagTwoSolve<LatticeFermionD> solver(mCG);

  /////////////////////////////////////////
  // Calculate Wall Source Propagators
  /////////////////////////////////////////


  // vector<Coordinate> point_srcs(1);  // FIXME: should iterate over many point sources 
  // point_srcs[0] = Coordinate(std::vector<int>{0,0,0,0});

  // int T = fdims[3];
  // vector<Coordinate> point_srcs(T);  
  // for(int i=0; i<T; ++i) point_srcs[i] = Coordinate(std::vector<int>{0,0,0,i});

  // vector<vector<int>> point_srcs = {{0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,0,3}, {0,0,0,4}, {0,0,0,5}, {0,0,0,6}, {0,0,0,7}};
  vector<vector<int>> point_srcs = {{1,0,0,0}, {1,0,0,1}, {1,0,0,2}, {1,0,0,3}};

  for(const auto &point_src: point_srcs) {
    std::cout << GridLogMessage << "point_src: " << point_src << std::endl;

    LatticePropagator prop(UGrid);
    for(int s = 0; s < 4; ++s) {
      for(int c = 0; c < 3; ++c) {
        std::cout << GridLogMessage << "s: " << s << " c: " << c << std::endl;
        LatticeFermion src(UGrid), src5d(FGrid), sol(UGrid), sol5d(FGrid);

        // Construct 4d source

        src = Zero();
        LatticeFermionSite site = Zero();
        site()(s)(c) = 1.;
        pokeSite(site, src, point_src);

        // Solve propagator
        Dmob.ImportPhysicalFermionSource(src, src5d);

        solver(Dmob, src5d, sol5d); // zyd: Dmob is not used, but required for syntax // In the MixedPrecisionConjugateGradientOp class, Dmob is not used

        Dmob.ExportPhysicalFermionSolution(sol5d, sol); // 5d->4d, v(x) = P_L v(x,0) + P_R v(x, Ls-1)

        FermToProp<WilsonImplR>(prop, sol, s, c);
      }
    }

    // string prop_fname = output_prefix + "/point_l/" + coor2str(point_src.toVector());
    // string prop_fname = output_prefix + "/point_" + quark_type + "/" + coor2str(point_src.toVector());
    string prop_fname = output_prefix + "/point_" + quark_type + "/" + coor2str(point_src);
    writeScidac_prop_d2f(prop, prop_fname);

    // // std::cout << prop << std::endl;
    // std::cout << "Full propagator at {0,0,0,0}" << std::endl;
    // print_grid_field_site(prop, {0,0,0,0});
    //
    //
    // LatticeSpinMatrix tmp = peekColour(prop, 0, 0); // The color index dependence of prop is delta_{a,b}, and is trivial
    // // std::cout << "Propagator without color index" << std::endl;
    // // print_grid_field_site(tmp, {0,1,2,3});
    // // print_grid_field_site(tmp, {3,2,0,1});
    // // print_grid_field_site(tmp, {3,2,1,0});
    //
    // std::cout << "Propagator without color index" << std::endl;
    // std::cout << tmp << std::endl;
  }  // end of loop of point src



  std::cout << "FINISHED!" << std::endl;
  Grid_finalize();
}

