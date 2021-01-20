// This file is the same as point_l.cc, except that it performs a random gauge transformation on Umu, and check if the resulting propagator following L(x, x0) -> g(x) L(x, x0) g(x0)^dagger
// 
// ??? How does ZMobius eigenvector change after gauge transformation??  In this code, I am not using eigenvectors. It takes about 0.5 h to do one CG, and about 6h for one point source.

#include "/direct/sdcc+u/ydzhao/A2AGrid/read_compressed.h"
#include "kaon/kaon.h" // do not know why, but kaon.h must be after a2a_field.h
#include "amplitude/form_factor.h"

// On 16 nodes, it takes ~10min to solve for one point source

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

void load_compressed_evecs(const std::string &evec_dir, std::vector<double> &evals, std::vector<LatticeFermionF> &evecs_f, Grid::GridCartesian *FGrid_f, Grid::GridRedBlackCartesian *FrbGrid_f) {

  USING_NAMESPACE_CPS

  if(evec_dir == "None") {
    std::cout << "Not using eigenvectors" << std::endl;
  }
  else {
    typedef A2ApoliciesDoubleAutoAlloc A2Apolicies;
    typedef A2Apolicies::FgridGFclass LatticeType;

    DoArg do_arg;
    if(!do_arg.Decode("do_arg.vml","do_arg")){
      do_arg.Encode("do_arg.templ","do_arg");
      VRB.Result("","main","Can't open do_arg.vml!\n");exit(1);
    }   

    int nthreads = omp_get_max_threads();
    std::cout << "Number of threads " << nthreads << std::endl;
    GJP.Initialize(do_arg);  // This sets the number of threads to 1
    GJP.SetNthreads(nthreads);
    std::cout << "Setting number of threads " << nthreads << std::endl;

    //Initialize FGrid
    FgridParams fgp;
    fgp.mobius_scale = 4.0; //b+c 
    LatticeType lattice(fgp);
    lattice.ImportGauge();

    std::cout << GridLogMessage << "Before reading compressed: " << std::endl;
    zyd_read_compressed(evec_dir, FGrid_f, FrbGrid_f, lattice, evals, evecs_f);

    std::cout << "Number of eigenvectors: " << evals.size()  << std::endl;
  }
}




void restrictTimeRange(LatticePropagator &lat, int vt)  { // The allowed interval of u is: vt <= ut <= vt+16

  const int T = lat.Grid()->_fdimensions[3];

  thread_for(ss, lat.Grid()->lSites(), {
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

    int ut = gcoor[3];

    int dist;   // dist = u_t - v_t, taking periodic bundary condition into account
    if(abs(ut - vt)<=T/2) dist = ut - vt;
    else if(abs(ut + T - vt) <= T/2) dist = ut + T - vt;
    else dist = ut - T - vt;

    LatticePropagatorSite m;
    if(dist==0) {}                      // if u_t == v_t , do nothing
    else if(dist > 0 && dist <= T/4) {  // u_t>v_t && u_t - v_t <= T/4  // if u is on the right of v, multiple it by 2
      peekLocalSite(m, lat, lcoor);
      m = 2. * m;
      pokeLocalSite(m, lat, lcoor);
    }   
    else {                           // else, set to 0
      m = Zero(); 
      pokeLocalSite(m, lat, lcoor);
    }   

  }); 
}




LatticeKGG calc_leptonic_with_coef(double M_K, const std::vector<int> &v, GridCartesian *grid) {  // return L_munu(u, v), where v is a fixed point.

  LatticeKGG lep(grid);

  form_factor_integrand(lep, M_K);

  double lep_coef = 2. / std::pow(M_K, 4);
  // double hadron_coef = env.Z_V * env.Z_V * 2. * env.M_K / env.N_K;  // Note: I did not multiply hadronic coefficient
  
  lep = lep * lep_coef;

  for(int mu=0; mu<4; ++mu) lep = Cshift(lep, mu, -v[mu]); // shift v to origin // FIXME: check it is right to shift -v[mu]

  return lep;
}








int main(int argc, char **argv)
{
  cps::Start(&argc, &argv); // Grid_init(&argc,&argv) is called inside this function

  int N_pts_sources = 50; // Number of point source to calculate in one job

  int Ls_outer = 24, Ls_inner = 12;
  double mass = 0.00107, b = 2.5, M5 = 1.8;  // !! mass must be the mass of light quark
  double M_K = 0.50365;
  vector<int> fdims = {24, 24, 24, 64};

  std::vector<std::complex<double>> omega;
  omega.resize(12);
  omega[0] = std::complex<double>(1.0903256131299373, 0); 
  omega[1] = std::complex<double>(0.9570283702230611, 0); 
  omega[2] = std::complex<double>(0.7048886040934104, 0); 
  omega[3] = std::complex<double>(0.48979921782791747, 0); 
  omega[4] = std::complex<double>(0.328608311201356, 0); 
  omega[5] = std::complex<double>(0.21664245377015995, 0); 
  omega[6] = std::complex<double>(0.14121112711957107, 0); 
  omega[7] = std::complex<double>(0.0907785101745156, 0); 
  omega[8] = std::complex<double>(0.05608303440064219, -0.007537158177840385);
  omega[9] = std::complex<double>(0.05608303440064219, 0.007537158177840385);
  omega[10] = std::complex<double>(0.0365221637144842, -0.03343945161367745);
  omega[11] = std::complex<double>(0.0365221637144842, 0.03343945161367745);

  // string output_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_my_props/sequential";
  string output_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_my_props/sequential_gauge_trans_test";
  string original_prop_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_my_props/sequential"; // propagator for Umu without gauge transformation
  string Umu_dir = "/hpcgpfs01/work/lqcd/etap/chulwoo/evec/24ID1Gev/configurations";
  string evec_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/evecs";
  // string point_l_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_my_props/point_l";
  string point_l_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_my_props/point_l_gauge_trans_test";

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

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(fdims, GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian *UGrid_f = SpaceTimeGrid::makeFourDimGrid(fdims, GridDefaultSimd(4, vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);

  // For inner ZMobius
  GridCartesian *FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls_inner, UGrid);
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls_inner, UGrid);
  GridCartesian *FGrid_f = SpaceTimeGrid::makeFiveDimGrid(Ls_inner, UGrid_f);
  GridRedBlackCartesian *FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls_inner, UGrid_f);

  // For outer Mobius
  GridCartesian *mobFGrid = SpaceTimeGrid::makeFiveDimGrid(Ls_outer, UGrid);
  GridRedBlackCartesian *mobFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls_outer, UGrid);
  GridCartesian *mobFGrid_f = SpaceTimeGrid::makeFiveDimGrid(Ls_outer, UGrid_f);
  GridRedBlackCartesian *mobFrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls_outer, UGrid_f);

  // genereate a random gauge transformation field
  LatticeColourMatrix g_test(UGrid);
  GridParallelRNG pRNG(UGrid); pRNG.SeedFixedIntegers({1,2,3,4});
  SU<3>::LieRandomize(pRNG, g_test, 1.0);  // g_test must be the same as the g_test in point_l_gauge_trans_test.cc

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {

    /////////////////////////////////////////
    // Find point sources to calculate
    /////////////////////////////////////////
    makeFileDir(output_prefix + "/" + to_string(traj) + "/xyz", UGrid); // remove file name, and make directory
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "makeFildDir " << output_prefix + "/" + to_string(traj) << std::endl;

    vector<vector<int>> new_pts;
    vector<vector<int>> all_pts = my_get_xgs(point_l_prefix + "/" + to_string(traj));
    vector<vector<int>> completed_pts = my_get_xgs(output_prefix + "/" + to_string(traj));
    for(const vector<int> &v: all_pts) {
      if(find(completed_pts.begin(), completed_pts.end(), v) == completed_pts.end()) new_pts.push_back(v);
      if(new_pts.size() >= N_pts_sources) break;
    }
    MPI_Barrier(MPI_COMM_WORLD);  // make sure every node has read complated_pts, before creating new files

    if(UGrid->IsBoss()) { // use an empty file to indicate that this propagator is being calculated by another job
      for(const vector<int> &v: new_pts) {
        ofstream f(output_prefix + "/" + to_string(traj) + "/" + coor2str(v)); 
        f.close();
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << GridLogMessage << "Will calculate " << new_pts.size() << " point sources in this job" << std::endl;

    // if(new_pts.empty()) continue; // For gauge transformation test, I commented this out.

    /////////////////////////////////////////
    // Read Gauge Field and do Gauge fixing
    /////////////////////////////////////////
    LatticeGaugeField Umu(UGrid);
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, Umu_dir + "/ckpoint_lat." + to_string(traj));

    // !!! perform gauge transformation
    SU<3>::GaugeTransform(Umu, g_test); // U_mu(x) -> g(x) U_mu(x) g(x+mu)^dagger

    /////////////////////////////////////////
    // Read eigenvectors 
    /////////////////////////////////////////
    string evec_dir = evec_prefix + "/" + to_string(traj) + "/lanczos.output";

    vector<double> evals;
    vector<LatticeFermionF> evecs_f;
    // load_compressed_evecs(evec_dir, evals, evecs_f, FGrid_f, FrbGrid_f);  // !!! MUST NOT use eigenvectors; the eigenvectors change after we do a gauge transformation on U

    /////////////////////////////////////////
    // Set Operators
    /////////////////////////////////////////

    LatticeGaugeFieldF Umu_f(UGrid_f);
    precisionChange(Umu_f, Umu);

    // Must set boundary phase in time direction to -1
    typename ZMobiusFermionD::ImplParams params;
    std::vector<Complex> boundary_phases(4, 1.);
    boundary_phases[3] = -1.;
    params.boundary_phases = boundary_phases;

    MobiusFermionD Dmob = MobiusFermionD(Umu, *mobFGrid, *mobFrbGrid, *UGrid, *UrbGrid, mass, M5, b, b-1., params);
    ZMobiusFermionF Dzmob_f = ZMobiusFermionF(Umu_f, *FGrid_f, *FrbGrid_f, *UGrid_f, *UrbGrid_f, mass, M5, omega, 1., 0., params);

    // CG_outer is used only for PV
    double resid = 1e-7, resid_outer = 1e-4, resid_inner = 1e-4;
    int max_iters = 10000;

    ConjugateGradient<LatticeFermionD> CG_outer(resid_outer, max_iters); 
    using PVtype = PauliVillarsSolverFourierAccel<LatticeFermionD, LatticeGaugeFieldD>;
    PVtype PV_outer(Umu, CG_outer);

    ConjugateGradient<LatticeFermionF> CG_inner(resid_inner, max_iters, 0); 
    using SchurSolverF = SchurRedBlackDiagTwoSolve<LatticeFermionF>;
    SchurSolverF SchurSolver_inner(CG_inner);

    using GuessF = DeflatedGuesser<LatticeFermionF>;
    GuessF guesser_inner(evecs_f, evals);

    MADWF<MobiusFermionD, ZMobiusFermionF, PVtype, SchurSolverF, GuessF> madwf(Dmob, Dzmob_f, PV_outer, SchurSolver_inner, guesser_inner, resid, 100);

    /////////////////////////////////////////
    // Calculate Sequential Propagators
    /////////////////////////////////////////
    
    new_pts.clear(); new_pts.push_back(std::vector<int>{0,0,22,7}); // Test Only One point source

    for(const vector<int> &v: new_pts) {
      cout << GridLogMessage << "v: " << v << endl;

      /////////////////////////////////////////
      // Construct Source
      /////////////////////////////////////////
      
      LatticeKGG lep = calc_leptonic_with_coef(M_K, v, UGrid);

      LatticePropagator Luv(UGrid);
      string point_l_path = point_l_prefix + "/" + to_string(traj) + "/" + coor2str(v);
      readScidac_prop_f2d(Luv, point_l_path);

      LatticePropagator fullSrc(UGrid);  fullSrc = Zero();
      for(int mu=0; mu<4; ++mu) {
        for(int nu=0; nu<4; ++nu) {
          if( mu==3 || nu==3 || mu==nu) continue;  // For these directions, lep_munu = 0
          LatticeComplex lep_munu = PeekIndex<LorentzIndex>(lep, mu, nu);
          fullSrc += gmu[mu] * Luv * gmu[nu] * lep_munu;
        }
      }

      restrictTimeRange(fullSrc, v[3]); // source is non-zero only for where vt <= ut <= vt+T/4; When ut>vt, multiply by 2

      /////////////////////////////////////////
      // Run CG (MADWF)
      /////////////////////////////////////////

      LatticePropagator prop(UGrid);
      for(int s = 0; s < 4; ++s) {
        for(int c = 0; c < 3; ++c) {
          std::cout << GridLogMessage << "s: " << s << " c: " << c << std::endl;

          LatticeFermion src(UGrid), sol(UGrid), sol5d(Dmob.FermionGrid());

          PropToFerm<WilsonImplR>(src, fullSrc, s, c);

          madwf(src, sol5d); // 4d src -> 5d solution

          Dmob.ExportPhysicalFermionSolution(sol5d, sol); // 5d->4d, v(x) = P_L v(x,0) + P_R v(x, Ls-1)

          FermToProp<WilsonImplR>(prop, sol, s, c);
        }
      }
      string prop_fname = output_prefix + "/" + to_string(traj)  + "/" + coor2str(v);
      writeScidac_prop_d2f(prop, prop_fname);

      // Test gauge transformation relationship: L(x, x0) -> g(x) L(x, x0) g(x0)^dagger 
      
      LatticePropagator prop2(UGrid);  // propagator calculated on the original lattice gauge
      string prop2_fname = original_prop_prefix + "/" + to_string(traj)  + "/" + coor2str(v);
      readScidac_prop_f2d(prop2, prop2_fname);

      LatticeColourMatrixSite g_test_x0;
      peekSite(g_test_x0, g_test, Coordinate(v));
      prop2 = g_test * prop2 * adj(g_test_x0); // Gauge transformation test: L(x,x0) -> g(x) L(x,x0) g(x0)^dagger

      print_grid_field_site(prop, {0,1,2,3});
      print_grid_field_site(prop2, {0,1,2,3});
      LatticePropagator tmp = prop - prop2;
      std::cout << "diff: " << norm2(tmp) << std::endl;  // Guage transformation test is passed if prop == prop2



    }  // end of loop of point sources

  }  // end of loop of traj


  std::cout << "Finished!" << std::endl;
  Grid_finalize();
}

