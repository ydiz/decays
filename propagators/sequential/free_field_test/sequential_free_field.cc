// This file is adapted from sequential.cc; However, the solvers are not MADWF, but the MixedPrecisionCG from propagator/point_s/point_s.cc

#include "kaon/kaon.h" // do not know why, but kaon.h must be after a2a_field.h

// On 16 nodes, it takes ~10min to solve for one point source

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



//
// LatticeKGG calc_leptonic_with_coef(double M_K, const std::vector<int> &v, GridCartesian *grid) {  // return L_munu(u, v), where v is a fixed point.
//
//   LatticeKGG lep(grid);
//
//   form_factor_integrand(lep, M_K);
//
//   double lep_coef = 2. / std::pow(M_K, 4);
//   // double hadron_coef = env.Z_V * env.Z_V * 2. * env.M_K / env.N_K;  // Note: I did not multiply hadronic coefficient
//
//   lep = lep * lep_coef;
//
//   for(int mu=0; mu<4; ++mu) lep = Cshift(lep, mu, -v[mu]); // shift v to origin 
//
//   return lep;
// }
//






int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  // int Ls_outer = 24, Ls_inner = 12;
  // double mass = 0.00107, b = 2.5, M5 = 1.8;  // !! mass must be the mass of light quark
  // double M_K = 0.50365;
  // vector<int> fdims = {24, 24, 24, 64};

  int Ls = 16;
  double mass = 0.04, b = 1.0, M5 = 1.0; // light quark
  vector<int> fdims = {8, 8, 8, 8};
  double M_K = 0.5;
  int max_uv_sep = 3; // The maximum separation between u and v
  // double time_boundary_phase = 0.; // No boundary phase; periodic boundary condition
  double time_boundary_phase = 113. * M_PI / 180.;  // 113 degrees

  std::cout << "Lat size: " << fdims << std::endl;
  std::cout << "Ls: " << Ls << std::endl;
  std::cout << "M5: " << M5 << std::endl;
  std::cout << "mass: " << mass << std::endl;
  std::cout << "b: " << b << std::endl;
  std::cout << "M_K: " << M_K << std::endl;
  std::cout << "max_uv_sep: " << max_uv_sep << std::endl;
  std::cout << "time_boundary_phase: " << time_boundary_phase << std::endl;


  string output_prefix = "./free_field_props";
  string point_l_prefix = "/home/ydzhao/cuth/decays/propagators/point_s/free_field_test/free_field_props/point_l"; // Free field propagator

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
  Umu = 1.0;
  // FieldMetaData header;
  // NerscIO::readConfiguration(Umu, header, Umu_dir + "/ckpoint_lat." + to_string(traj));

  /////////////////////////////////////////
  // Set Operators
  /////////////////////////////////////////

  LatticeGaugeFieldF Umu_f(UGrid_f);
  precisionChange(Umu_f, Umu);

  typename MobiusFermionD::ImplParams params;
  std::vector<Complex> boundary_phases(4, 1.);  // For free field test, all directions must be periodic
  boundary_phases[3] = std::exp(Complex(0., time_boundary_phase)); // exp(i alpha)
  // boundary_phases[3] = -1.;   // 
  params.boundary_phases = boundary_phases;



  MobiusFermionD Dmob = MobiusFermionD(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, b, b-1., params);
  MobiusFermionF Dmob_f = MobiusFermionF(Umu_f, *FGrid_f, *FrbGrid_f, *UGrid_f, *UrbGrid_f, mass, M5, b, b-1., params);
  SchurDiagTwoOperator<MobiusFermionD, LatticeFermionD> mobHermOp = SchurDiagTwoOperator<MobiusFermionD, LatticeFermionD>(Dmob);
  SchurDiagTwoOperator<MobiusFermionF, LatticeFermionF> mobHermOp_f = SchurDiagTwoOperator<MobiusFermionF, LatticeFermionF>(Dmob_f);

  // double resid = 1e-8;
  double resid = 1e-16;
  int max_iters = 10000;
  MixedPrecisionConjugateGradientOp<LatticeFermionD, LatticeFermionF> mCG(resid, max_iters, 50, FrbGrid_f, mobHermOp_f, mobHermOp);
  SchurRedBlackDiagTwoSolve<LatticeFermionD> solver(mCG);

  /////////////////////////////////////////
  // Calculate Sequential Propagators
  /////////////////////////////////////////

  // vector<vector<int>> new_pts = {{0,0,0,0}};
  vector<vector<int>> new_pts = {{0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,0,3}, {0,0,0,4}, {0,0,0,5}, {0,0,0,6}, {0,0,0,7}};
  for(const vector<int> &v: new_pts) {
    cout << GridLogMessage << "v: " << v << endl;

    /////////////////////////////////////////
    // Construct Source
    /////////////////////////////////////////

    // LatticeKGG lep = calc_leptonic_with_coef(M_K, v, UGrid);
    LatticeKGG lep(UGrid);
    EM_factor_half_lattice(lep, v, M_K, max_uv_sep);

    LatticePropagator Luv(UGrid);
    string point_l_path = point_l_prefix + "/" + coor2str(v);
    readScidac_prop_f2d(Luv, point_l_path);

    // std::cout << "Luv:" << std::endl;
    // std::cout << Luv << std::endl;
    // return 0;

    LatticePropagator fullSrc(UGrid);  fullSrc = Zero();
    for(int mu=0; mu<4; ++mu) {
      for(int nu=0; nu<4; ++nu) {
        if( mu==3 || nu==3 || mu==nu) continue;  // For these directions, lep_munu = 0
        LatticeComplex lep_munu = PeekIndex<LorentzIndex>(lep, mu, nu);
        fullSrc += gmu[mu] * Luv * gmu[nu] * lep_munu;
      }
    }

    // std::cout << "fullSrc:" << std::endl;
    // std::cout << fullSrc << std::endl;
    // return 0;

    // restrictTimeRange(fullSrc, v[3]); // source is non-zero only for where vt <= ut <= vt+T/4; When ut>vt, multiply by 2

    /////////////////////////////////////////
    // Run CG (MADWF)
    /////////////////////////////////////////

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

        // LatticeFermion src(UGrid), sol(UGrid), sol5d(Dmob.FermionGrid());
        //
        // PropToFerm<WilsonImplR>(src, fullSrc, s, c);
        //
        // madwf(src, sol5d); // 4d src -> 5d solution
        //
        // Dmob.ExportPhysicalFermionSolution(sol5d, sol); // 5d->4d, v(x) = P_L v(x,0) + P_R v(x, Ls-1)
        //
        // FermToProp<WilsonImplR>(prop, sol, s, c);
      }
    }
    string prop_fname = output_prefix + "/" + coor2str(v);
    writeScidac_prop_d2f(prop, prop_fname);

    // // std::cout << prop << std::endl;
    // std::cout << "Full propagator at {0,0,0,0}" << std::endl;
    // print_grid_field_site(prop, {0,0,0,0});
    //
    //
    // LatticeSpinMatrix tmp = peekColour(prop, 0, 0); // The color index dependence of prop is delta_{a,b}, and is trivial
    // std::cout << "Propagator without color index" << std::endl;
    // std::cout << tmp << std::endl;

  }  // end of loop of point sources



  std::cout << "Finished!" << std::endl;
  Grid_finalize();
}

