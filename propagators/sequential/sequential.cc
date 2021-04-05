// !! Old version Grid has a bug. In Grid/qcd/action/fermion/FourierAcceleratedPV.h, the boundary phase of WilsonTMFermion5D is not set. 
// !! If this bug is not corrected, MADWF does not converge.

#include "kaon/kaon.h" 
#include "../read_compressed_evec.h"

// At BNL: On 16 nodes, it takes ~10min to solve for one point source
// At Cori: On 32 nodes, it takes ~800s to solve for one point source, and ~11h for 50 point sources

using namespace std;
using namespace Grid;

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);
  // cps::Start(&argc, &argv); // Grid_init(&argc,&argv) is called inside this function

  // int N_pts_sources = 1;  // FIXME: for test
  // int N_pts_sources = 50; // Number of point source to calculate in one job
  int N_pts_sources = 100; // Number of point source to calculate in one job

  string output_prefix = "/global/cfs/cdirs/mp13/ydzhao/24ID/sequential";
  // string output_prefix = "."; // for test
  string Umu_dir = "/global/cfs/cdirs/mp13/ydzhao/24ID/configurations";
  string evec_prefix = "/global/cscratch1/sd/ydzhao/evecs";
  string point_l_prefix = "/global/cscratch1/sd/ydzhao/point_l";

  int Ls_outer = 24, Ls_inner = 12;
  double mass = 0.00107, b = 2.5, M5 = 1.8;  // !! mass must be the mass of light quark
  double M_K = 0.50365;
  vector<int> fdims = {24, 24, 24, 64};


  int max_uv_sep = 16;   // The maximum separation between u and v to sum over


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

    if(new_pts.empty()) continue;

    /////////////////////////////////////////
    // Read Gauge Field and do Gauge fixing
    /////////////////////////////////////////
    LatticeGaugeField Umu(UGrid);
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, Umu_dir + "/ckpoint_lat." + to_string(traj));

    /////////////////////////////////////////
    // Read eigenvectors 
    /////////////////////////////////////////
    string evec_dir = evec_prefix + "/" + to_string(traj) + "/lanczos.output";

    vector<double> evals;
    vector<LatticeFermionF> evecs_f;
    int ngroups = 1;  // To save memory, Should be number of MPI processes per node
    zyd_read_compressed_evecs(evec_dir, FrbGrid_f, evecs_f, evals, ngroups); 
    // // load_compressed_evecs(evec_dir, evals, evecs_f, FGrid_f, FrbGrid_f);  

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
    
    // for(const vector<int> &v: new_pts) {
    for(int i=0; i<new_pts.size(); ++i) {

      vector<int> v = new_pts[i];

      cout << GridLogMessage << "Number of points: " << i << endl;
      cout << GridLogMessage << "v: " << v << endl;

      /////////////////////////////////////////
      // Construct Source
      /////////////////////////////////////////
      
      // LatticeKGG lep = calc_leptonic_with_coef(M_K, v, UGrid);
      LatticeKGG lep(UGrid);
      EM_factor_half_lattice(lep, v, M_K, max_uv_sep);

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

      // restrictTimeRange(fullSrc, v[3]); // source is non-zero only for where vt <= ut <= vt+T/4; When ut>vt, multiply by 2

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
    }  // end of loop of point sources

  }  // end of loop of traj


  std::cout << "Finished!" << std::endl;
  Grid_finalize();
}

