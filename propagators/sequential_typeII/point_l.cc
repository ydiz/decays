#include "/direct/sdcc+u/ydzhao/A2AGrid/read_compressed.h"
#include "../../kaon/kaon.h" // do not know why, but kaon.h must be after a2a_field.h

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
    fgp.mobius_scale = 4.0; //b+c //FIXME: read from argument
    LatticeType lattice(fgp);
    lattice.ImportGauge();

    std::cout << GridLogMessage << "Before reading compressed: " << std::endl;
    zyd_read_compressed(evec_dir, FGrid_f, FrbGrid_f, lattice, evals, evecs_f);

    std::cout << "Number of eigenvectors: " << evals.size()  << std::endl;
  }
}









int main(int argc, char **argv)
{
  cps::Start(&argc, &argv); // Grid_init(&argc,&argv) is called inside this function

  int Ls_outer = 24, Ls_inner = 12;
  double mass = 0.00107, b = 2.5, M5 = 1.8;  // !! mass must be the mass of light quark
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

  string output_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_my_props";
  string Umu_dir = "/hpcgpfs01/work/lqcd/etap/chulwoo/evec/24ID1Gev/configurations";
  string evec_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/evecs";

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
    load_compressed_evecs(evec_dir, evals, evecs_f, FGrid_f, FrbGrid_f);

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
    // Calculate Wall Source Propagators
    /////////////////////////////////////////

    vector<Coordinate> point_srcs(1);  // FIXME: should iterate over many point sources 
    point_srcs[0] = Coordinate(std::vector<int>{0,0,7,25});

    for(Coordinate point_src: point_srcs) {
      std::cout << GridLogMessage << "point_src: " << point_src << std::endl;

      LatticePropagator prop(UGrid);
      for(int s = 0; s < 4; ++s) {
        for(int c = 0; c < 3; ++c) {
          std::cout << GridLogMessage << "s: " << s << " c: " << c << std::endl;

          LatticeFermion src(UGrid), sol(UGrid), sol5d(Dmob.FermionGrid());

          src = Zero();
          LatticeFermionSite site = Zero();
          site()(s)(c) = 1.;
          pokeSite(site, src, point_src);

          madwf(src, sol5d); // 4d src -> 5d solution

          Dmob.ExportPhysicalFermionSolution(sol5d, sol); // 5d->4d, v(x) = P_L v(x,0) + P_R v(x, Ls-1)

          FermToProp<WilsonImplR>(prop, sol, s, c);
        }
      }

      string prop_fname = output_prefix + "/my_point_l_test/" + to_string(traj)  + "/" + coor2str(point_src.toVector());
      writeScidac_prop_d2f(prop, prop_fname);
    }  // end of loop of tW

  }  // end of loop of traj


  Grid_finalize();
}

