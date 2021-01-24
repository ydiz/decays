#include "kaon/kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


vector<LatticePropagator> conv_with_E(const vector<LatticePropagator> &F, double M_K, int max_uv_sep) {
  using namespace std;
  GridBase *grid = F[0].Grid();

  FFT theFFT((GridCartesian *)grid);

  // static vector<vector<LatticeComplex>> Euv_fft_conj(4, vector<LatticeComplex>(4, grid));  // // Euv_fft = (Fourier{Euv})^*, where Euv is EM factor 
  static vector<vector<LatticeComplex>> Euv_fft_conj;  // Euv_fft = (Fourier{Euv})^*, where Euv is EM factor 
  static bool initialized = false;
  if(!initialized) {
    Euv_fft_conj = calc_Euv_fft_conj(grid, M_K, max_uv_sep);
    initialized = true;
  }

  vector<LatticePropagator> F_fft(3, grid); // F_fft[3] is not used; so I do not calculate it
  for(int mu=0; mu<3; ++mu) theFFT.FFT_all_dim(F_fft[mu], F[mu], FFT::forward);

  vector<LatticePropagator> rst(4, grid);
  for(int mu=0; mu<4; ++mu) rst[mu] = Zero();
  for(int nu=0; nu<3; ++nu) {   // Here, nu != 3, because E_{3, mu} = 0
    for(int mu=0; mu<3; ++mu) { // Here, mu != 3, because E_{nu, 3} = 0
      if(mu==nu) continue;  // if mu==nu, E_{mu,mu} = 0
      rst[nu] += Euv_fft_conj[nu][mu] * F_fft[mu];
    }
  }

  for(int nu=0; nu<3; ++nu) theFFT.FFT_all_dim(rst[nu], rst[nu], FFT::backward); // rst[3] is always 0

  return rst;
}

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);
  // zyd_init_Grid_Qlattice(argc, argv);
 
 
#ifdef CUTH_FREE_FIELD
  int tsep = 1;
  int max_uv_sep = 3;
  Env env("FreeField_8nt8");
#else
  int tsep = 16;
  int max_uv_sep = 16;
  Env env("24ID");
#endif

  // env.N_pt_src = 1;  // keep only one point
  env.N_pt_src = -1;  // Use all points

  /////////////////////////////////////////////////

  const int T = env.grid->_fdimensions[3];

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

  // int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    vector<LatticePropagator> wl = env.get_wall('l');
    vector<LatticePropagator> ws = env.get_wall('s');

    vector<Complex> rst_Q1_vt_avgSrc(T, 0.), rst_Q2_vt_avgSrc(T, 0.);

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_s.resize(env.N_pt_src);
    for(const auto &x: env.xgs_s) {
      std::cout << "# Point source: " << num_pt_src << std::endl;
      std::cout << "Point source: " << x << std::endl;
      ++num_pt_src;

      LatticePropagator pl = env.get_point(x, 'l'); // pl = L(x, v) or L(u, v)
      LatticePropagator ps = env.get_point(x, 's'); // ps = H(x, v) or H(u, v)
      LatticePropagatorSite Lxx; // L(x, x)
      peekSite(Lxx, pl, x);

      int tK = (x[3] - tsep + T) % T;

      vector<LatticePropagator> Fu(4, env.grid), Gv(4, env.grid); // F_mu(u, x)
      for(int mu=0; mu<4; ++mu) Fu[mu] = adj(pl) * gmu5[mu] * wl[tK];

      for(int nu=0; nu<4; ++nu) Gv[nu] = adj(ws[tK]) * gmu5[nu] * ps; 
      LatticeComplex exp_factor = exp_v0_tK(env.grid, tK, env.M_K);
      for(int nu=0; nu<4; ++nu) Gv[nu] = Gv[nu] * exp_factor;           // G_nu(v, x) *= exp(M_k * (v_0 - t_K))

      vector<LatticePropagator> Cv = conv_with_E(Fu, env.M_K, max_uv_sep);

      LatticePropagator A(env.grid); A = Zero();
      for(int nu=0; nu<4; ++nu) A += Cv[nu] * Gv[nu];
      A = A - adj(A);  // Contribution of K0_bar
      // std::cout << A << std::endl;
      // exit(0);

      LatticeComplex rst_Q1(env.grid); rst_Q1 = Zero();
      LatticeComplex rst_Q2(env.grid); rst_Q2 = Zero();
      for(int rho=0; rho<4; ++rho) {
        rst_Q1 += trace(gL[rho] * Lxx) * trace(gL[rho] * A);
        rst_Q2 += trace(gL[rho] * Lxx * gL[rho] * A);
      }

      vector<LatticeComplexSite> rst_Q1_vt, rst_Q2_vt; // Sum over each time slice of v
      sliceSum(rst_Q1, rst_Q1_vt, Tdir);
      sliceSum(rst_Q2, rst_Q2_vt, Tdir);
      for(int i=0; i<T; ++i) {
        rst_Q1_vt_avgSrc[i] += rst_Q1_vt[(i+x[3])%T]()()(); // rst_Q1_avgSrc[0] is the time slice of point v
        rst_Q2_vt_avgSrc[i] += rst_Q2_vt[(i+x[3])%T]()()();
      }

    } // end of point source loop

    std::cout << GridLogMessage << "Number of point sources: " << num_pt_src << std::endl;
    for(int i=0; i<T; ++i) {
      rst_Q1_vt_avgSrc[i] *= (1. / double(num_pt_src)); // rst_Q1_avgSrc[0] is the time slice of point v
      rst_Q2_vt_avgSrc[i] *= (1. / double(num_pt_src));
    }

    std::cout << "traj [" << traj << "] amplitude Q1 by time slice: " << rst_Q1_vt_avgSrc << std::endl;
    std::cout << "traj [" << traj << "] amplitude Q2 by time slice: " << rst_Q2_vt_avgSrc << std::endl;

  } // end of traj loop

  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}
