#include "../../kaon.h"
#include "../control_R.h"

using namespace std;
using namespace Grid;

#ifdef CUTH_FREE_FIELD
#define NUM_XT 4
#define NUM_R 4
#else
#define NUM_XT 19
#define NUM_R 25
#endif


// return a vector of shape (N_t_seps, T, max_R) -> [t_sep][x_t][R]
vector<vector<vector<Complex>>> init_3d_vec(int N_t_seps, const Coordinate &fdims) {

  // const int T = fdims[3];
  // int X = fdims[0], Y = fdims[1], Z = fdims[2];
  // int num_R = 1 + int(calc_3d_vec_norm(X/2, Y/2, Z/2));
  vector<vector<vector<Complex>>> diagram(N_t_seps);
  for(auto &x: diagram) {
    // x.resize(T);
    x.resize(NUM_XT);
    for(auto &y: x) {
      // y.resize(num_R, 0.);
      y.resize(NUM_R, 0.);
    }
  }
  return diagram;
}

int main(int argc, char* argv[])
{
  // vector<Complex> a, b;
  // a+b;

  Grid_init(&argc, &argv);
  // zyd_init_Grid_Qlattice(argc, argv);

#ifdef CUTH_FREE_FIELD
  vector<int> t_seps = {1,2,3}; 
  Env env("FreeField_8nt8");
#else
  vector<int> t_seps = {12, 14, 16, 18};  // For eta: 10-18; for pion 18-24  // separation of eta and kaon walls
  Env env("24ID");
#endif

  const int T = env.grid->_fdimensions[3];

  int N_t_seps = t_seps.size();
  cout << "t_seps: " << t_seps << endl;

  int target_traj;
  if( GridCmdOptionExists(argv, argv+argc, "--traj") ) {
    string arg = GridCmdOptionPayload(argv, argv+argc, "--traj");
    GridCmdOptionInt(arg, target_traj);
  }
  else {
    cout << "--traj not specified; exiting" << endl;
    exit(0);
  }

  int traj_start = target_traj;
  int traj_end = target_traj;
  int traj_sep = 10;

  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  cout << string(20, '*') << endl;
  cout << "traj_start: " << traj_start << endl;
  cout << "traj_end: " << traj_end << endl;
  cout << "traj_sep: " << traj_sep << endl;
  cout << "traj_num: " << traj_num << endl;
  cout << string(20, '*') << endl;

  vector<string> diagrams = {"sBar_d_T1D1", "sBar_d_T1D2", "sBar_d_T2D1", "sBar_d_T2D2", "Hw_T1D1Q1", "Hw_T1D1Q2", "Hw_T2D1Q1", "Hw_T2D1Q2", "Hw_T2D2Q1", "Hw_T2D2Q2", "Hw_T3D1Q1", "Hw_T3D1Q2", "Hw_T3D2Q1", "Hw_T3D2Q2"};

  // map<string, decltype(&sBar_d_T1D1)> diagram_funcs {{"sBar_d_T1D1", sBar_d_T1D1}, {"sBar_d_T1D2", sBar_d_T1D2}, {"sBar_d_T2D1", sBar_d_T2D1}, {"sBar_d_T2D2", sBar_d_T2D2}, {"Hw_T1D1Q1", Hw_T1D1Q1}, {"Hw_T1D1Q2", Hw_T1D1Q2}, {"Hw_T2D1Q1", Hw_T2D1Q1}, {"Hw_T2D1Q2", Hw_T2D1Q2}, {"Hw_T2D2Q1", Hw_T2D2Q1}, {"Hw_T2D2Q2", Hw_T2D2Q2}, {"Hw_T3D1Q1", Hw_T3D1Q1}, {"Hw_T3D1Q2", Hw_T3D1Q2}, {"Hw_T3D2Q1", Hw_T3D2Q1}, {"Hw_T3D2Q2", Hw_T3D2Q2}};


  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    map<string, vector<vector<vector<Complex>>>> table_3d; 
    for(const string &diagram: diagrams) table_3d[diagram] = init_3d_vec(N_t_seps, env.grid->_fdimensions);

    // Read gauge transformation matrices for this trajectory
    LatticePropagator Lxx = env.get_Lxx();

    vector<LatticePropagator> wl = env.get_wall('l');
    vector<LatticePropagator> ws = env.get_wall('s');


    int num_pt_src = 0;
    for(const auto &v: env.xgs_l) {  // v is position of eta
      cout << "# Point source: " << num_pt_src << endl;
      cout << "Point source: " << v << endl;
      ++num_pt_src;

      LatticePropagator pl_eta = env.get_point(v, 'l'); // l = L(x, x_eta)
      LatticePropagator ps_eta = env.get_point(v, 's'); // l = L(x, x_eta)

      LatticePropagatorSite pl_eta_eta; peekSite(pl_eta_eta, pl_eta, v); 
      LatticePropagatorSite ps_eta_eta; peekSite(ps_eta_eta, ps_eta, v); 

      int t_v = v[3];
      for(int t_sep_idx=0; t_sep_idx<t_seps.size(); ++t_sep_idx) { // iterate through all possible separation between eta and kaon
        int t_sep = t_seps[t_sep_idx];
        int t_K = (t_v - t_sep + T) % T;

        // cout << GridLogMessage << "tK: " << t_K << endl;
        const LatticePropagator &wl_K = wl[t_K]; // L(x, t_K)
        const LatticePropagator &ws_K = ws[t_K]; // H(x, t_K)

        LatticePropagatorSite wl_eta_tK; peekSite(wl_eta_tK, wl_K, v); 
        LatticePropagatorSite ws_eta_tK; peekSite(ws_eta_tK, ws_K, v); 

        // sBar_d_T1D1 and Hw_T2D1
        LatticePropagator tmp = wl_eta_tK * adj(ws_K) * g5 * pl_eta;
        LatticeComplex sBar_d_T1D1 = trace(g5 * tmp);
        sBar_d_T1D1 = 2. * real(sBar_d_T1D1); // Add contribution of K0 bar

        tmp = tmp + adj(tmp); // Add contribution of K0 bar
        LatticeComplex Hw_T2D1Q1(pl_eta.Grid()); Hw_T2D1Q1 = Zero();
        for(int mu=0; mu<4; ++mu) Hw_T2D1Q1 += trace(gL[mu] * Lxx) * trace(gL[mu] * tmp);
        LatticeComplex Hw_T2D1Q2(pl_eta.Grid()); Hw_T2D1Q2 = Zero();
        for(int mu=0; mu<4; ++mu) Hw_T2D1Q2 += trace(gL[mu] * Lxx * gL[mu] * tmp);

        // sBar_d_T1D2 and Hw_T2D2
        tmp = wl_K * adj(ws_eta_tK) * g5 * adj(ps_eta);
        LatticeComplex sBar_d_T1D2 = trace(g5 * tmp);
        sBar_d_T1D2 = 2. * real(sBar_d_T1D2); // Add contribution of K0 bar

        tmp = tmp + adj(tmp); // Add contribution of K0 bar
        LatticeComplex Hw_T2D2Q1(pl_eta.Grid()); Hw_T2D2Q1 = Zero();
        for(int mu=0; mu<4; ++mu) Hw_T2D2Q1 += trace(gL[mu] * Lxx) * trace(gL[mu] * tmp);
        LatticeComplex Hw_T2D2Q2(pl_eta.Grid()); Hw_T2D2Q2 = Zero();
        for(int mu=0; mu<4; ++mu) Hw_T2D2Q2 += trace(gL[mu] * Lxx * gL[mu] * tmp);

        // sBar_d_T2D1/sBar_d_T2D2 and Hw_T3D1/Hw_T3D2
        tmp = wl_K * adj(ws_K);
        LatticeComplex sBar_d_T2D1 = trace(g5 * pl_eta_eta) * ( 2. * real(trace(g5 * tmp)) );
        LatticeComplex sBar_d_T2D2 = trace(g5 * ps_eta_eta) * ( 2. * real(trace(g5 * tmp)) );

        tmp = tmp + adj(tmp); // Add contribution of K0 bar
        LatticeComplex Hw_T3D1Q1(pl_eta.Grid()); Hw_T3D1Q1 = Zero();
        for(int mu=0; mu<4; ++mu) Hw_T3D1Q1 += trace(g5 * pl_eta_eta ) * trace(gL[mu] * Lxx) * trace(gL[mu] * tmp);
        LatticeComplex Hw_T3D1Q2(pl_eta.Grid()); Hw_T3D1Q2 = Zero();
        for(int mu=0; mu<4; ++mu) Hw_T3D1Q2 += trace(g5 * pl_eta_eta ) * trace(gL[mu] * Lxx * gL[mu] * tmp);
        LatticeComplex Hw_T3D2Q1(pl_eta.Grid()); Hw_T3D2Q1 = Zero();
        for(int mu=0; mu<4; ++mu) Hw_T3D2Q1 += trace(g5 * ps_eta_eta ) * trace(gL[mu] * Lxx) * trace(gL[mu] * tmp);
        LatticeComplex Hw_T3D2Q2(pl_eta.Grid()); Hw_T3D2Q2 = Zero();
        for(int mu=0; mu<4; ++mu) Hw_T3D2Q2 += trace(g5 * ps_eta_eta ) * trace(gL[mu] * Lxx * gL[mu] * tmp);

        // Hw_T1D1
        LatticePropagator tmp_eta = pl_eta * adj(pl_eta);
        LatticePropagator tmp_K = wl_K * adj(ws_K);
        tmp_K = tmp_K + adj(tmp_K); // Add contribution of K0 bar 

        LatticeComplex Hw_T1D1Q1(pl_eta.Grid()); Hw_T1D1Q1 = Zero();
        for(int mu=0; mu<4; ++mu) Hw_T1D1Q1 += trace(gL[mu] * tmp_eta) * trace(gL[mu] * tmp_K);
        LatticeComplex Hw_T1D1Q2(pl_eta.Grid()); Hw_T1D1Q2 = Zero();
        for(int mu=0; mu<4; ++mu) Hw_T1D1Q2 += trace(gL[mu] * tmp_eta * gL[mu] * tmp_K);

        map<string, LatticeComplex*> m = {{"sBar_d_T1D1", &sBar_d_T1D1}, {"sBar_d_T1D2", &sBar_d_T1D2}, {"sBar_d_T2D1", &sBar_d_T2D1}, {"sBar_d_T2D2", &sBar_d_T2D2}, {"Hw_T1D1Q1", &Hw_T1D1Q1}, {"Hw_T1D1Q2", &Hw_T1D1Q2}, {"Hw_T2D1Q1", &Hw_T2D1Q1}, {"Hw_T2D1Q2", &Hw_T2D1Q2}, {"Hw_T2D2Q1", &Hw_T2D2Q1}, {"Hw_T2D2Q2", &Hw_T2D2Q2}, {"Hw_T3D1Q1", &Hw_T3D1Q1}, {"Hw_T3D1Q2", &Hw_T3D1Q2}, {"Hw_T3D2Q1", &Hw_T3D2Q1}, {"Hw_T3D2Q2", &Hw_T3D2Q2}};

        for(const string &diagram: diagrams) {
          vector<vector<Complex>> rst = sum_T_R(*(m[diagram]), v); // rst[xt][R] // shape (T, max_R)
          // for(int xt=0; xt<rst.size(); ++xt) 
          //   for(int R=0; R<rst[0].size(); ++R) {
          for(int xt=0; xt<NUM_XT; ++xt) 
            for(int R=0; R<NUM_R; ++R) {
              table_3d[diagram][t_sep_idx][xt][R] += rst[(t_K+xt)%T][R]; // shift tK to 0
            }
        }

      }
    }

    for(const string &diagram: diagrams) {
      
      int num_xt = table_3d[diagram][0].size();
      int max_R = table_3d[diagram][0][0].size();
      for(int sep=0; sep<N_t_seps; ++sep) // average over number of point sources
        for(int xt=0; xt<num_xt; ++xt) 
          for(int R=0; R<max_R; ++R) 
            table_3d[diagram][sep][xt][R] /= double(num_pt_src);

      cout << "traj [" << traj << "] table3d " << diagram << ": " << table_3d[diagram] << endl; // table3d[t_sep][xt][R]
    }
  }

  cout << GridLogMessage << "Finished!" << endl;

  Grid_finalize();


  return 0;
}
