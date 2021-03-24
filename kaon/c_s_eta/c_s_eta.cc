#include "../kaon.h"
#include "c_s_eta_diagrams.h"

using namespace std;
using namespace Grid;

// return a vector of shape (traj_num, N_t_seps, T, T) -> [traj][t_sep][t_K][x_t]
std::vector<std::vector<std::vector<std::vector<Grid::Complex>>>> init_4d_vec(int traj_num, int N_t_seps, int T) {
  std::vector<std::vector<std::vector<std::vector<Grid::Complex>>>> diagram;
  diagram.resize(traj_num);
  for(auto &x: diagram) {
    x.resize(N_t_seps);
    for(auto &y: x) {
      y.resize(T);
      for(auto &z: y) z.resize(T);
    }
  }
  return diagram;
}


std::vector<std::vector<std::vector<Grid::Complex>>> average_tK(const std::vector<std::vector<std::vector<std::vector<Grid::Complex>>>> &diagram) {
  int traj_num = diagram.size();
  int sep_num = diagram[0].size();
  int tK_num = diagram[0][0].size();
  int xt_num = diagram[0][0][0].size();
  std::vector<std::vector<std::vector<Grid::Complex>>> rst(traj_num);
  for(int traj=0; traj<traj_num; ++traj) {
    rst[traj].resize(sep_num);
    for(int sep=0; sep<sep_num; ++sep) {
      rst[traj][sep].resize(xt_num);
      for(int xt=0; xt<xt_num; ++xt) {
        for(int tK=0; tK<tK_num; ++tK) {
          rst[traj][sep][xt] += diagram[traj][sep][tK][xt];
        }
        rst[traj][sep][xt] /= double(tK_num);
      }
    }
  }
  return rst;
}

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  // zyd_init_Grid_Qlattice(argc, argv);

#ifdef CUTH_FREE_FIELD
  vector<int> t_seps = {1,2,3}; 
  Env env("FreeField_8nt8");
#else
  vector<int> t_seps = {10, 12, 14, 16, 18, 20, 22, 24};  // For eta: 10-18; for pion 18-24  // separation of eta and kaon walls
  Env env("24ID");
#endif

  const int T = env.grid->_fdimensions[3];

  int N_t_seps = t_seps.size();
  std::cout << "t_seps: " << t_seps << std::endl;

  int target_traj;
  if( GridCmdOptionExists(argv, argv+argc, "--traj") ) {
    string arg = GridCmdOptionPayload(argv, argv+argc, "--traj");
    GridCmdOptionInt(arg, target_traj);
  }
  else {
    std::cout << "--traj not specified; exiting" << std::endl;
    exit(0);
  }

  int traj_start = target_traj;
  int traj_end = target_traj;
  int traj_sep = 10;

  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  vector<string> diagrams = {"sBar_d_T1D1", "sBar_d_T1D2", "sBar_d_T2D1", "sBar_d_T2D2", "Hw_T1D1Q1", "Hw_T1D1Q2", "Hw_T2D1Q1", "Hw_T2D1Q2", "Hw_T2D2Q1", "Hw_T2D2Q2", "Hw_T3D1Q1", "Hw_T3D1Q2", "Hw_T3D2Q1", "Hw_T3D2Q2"};

  map<string, decltype(&sBar_d_T1D1)> diagram_funcs {{"sBar_d_T1D1", sBar_d_T1D1}, {"sBar_d_T1D2", sBar_d_T1D2}, {"sBar_d_T2D1", sBar_d_T2D1}, {"sBar_d_T2D2", sBar_d_T2D2}, {"Hw_T1D1Q1", Hw_T1D1Q1}, {"Hw_T1D1Q2", Hw_T1D1Q2}, {"Hw_T2D1Q1", Hw_T2D1Q1}, {"Hw_T2D1Q2", Hw_T2D1Q2}, {"Hw_T2D2Q1", Hw_T2D2Q1}, {"Hw_T2D2Q2", Hw_T2D2Q2}, {"Hw_T3D1Q1", Hw_T3D1Q1}, {"Hw_T3D1Q2", Hw_T3D1Q2}, {"Hw_T3D2Q1", Hw_T3D2Q1}, {"Hw_T3D2Q2", Hw_T3D2Q2}};

  map<string, vector<vector<vector<vector<Complex>>>>> table_4d; 
  for(const string &diagram: diagrams) table_4d[diagram] = init_4d_vec(traj_num, N_t_seps, T);

  int traj_idx = 0;
  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    // Read gauge transformation matrices for this trajectory
    LatticeColourMatrix gt = env.get_gaugeTransform();
    LatticePropagator Lxx = env.get_Lxx();

    vector<LatticePropagator> wl = env.get_wall('l');
    vector<LatticePropagator> ws = env.get_wall('s');


    for(int t_K=0; t_K<T; ++t_K) { // iterate through the position of Kaon wall
      std::cout << GridLogMessage << "tK: " << t_K << std::endl;
      const LatticePropagator &wl_K = wl[t_K]; // L(x, t_K)
      const LatticePropagator &ws_K = ws[t_K]; // H(x, t_K)

      // Calculate L(t_eta, t_K); the sink must be in Coulomb gauge
      std::vector<LatticePropagatorSite> wl_K_sliceSum;
      {
        LatticePropagator wl_K_Coulomb = env.toCoulombSink(gt, wl_K);
        sliceSum(wl_K_Coulomb, wl_K_sliceSum, Tdir);
      }

      // Calculate H(t_eta, t_K); the sink must be in Coulomb gauge
      std::vector<LatticePropagatorSite> ws_K_sliceSum;
      {
        LatticePropagator ws_K_Coulomb = env.toCoulombSink(gt, ws_K);
        sliceSum(ws_K_Coulomb, ws_K_sliceSum, Tdir);
      }

      for(int t_sep_idx=0; t_sep_idx<t_seps.size(); ++t_sep_idx) { // iterate through all possible separation between eta and kaon
        int t_sep = t_seps[t_sep_idx];
        int t_eta = (t_K + t_sep) % T;

        const LatticePropagator &wl_eta = wl[t_eta]; // L(x, t_eta)
        const LatticePropagator &ws_eta = ws[t_eta]; // H(x, t_eta)
        LatticePropagatorSite wl_tEta_tK = wl_K_sliceSum[t_eta];  // L(t_eta, t_K)
        LatticePropagatorSite ws_tEta_tK = ws_K_sliceSum[t_eta];  // H(t_eta, t_K)
        LatticePropagatorSite wl_tEta_tEta, ws_tEta_tEta;
        {
          std::vector<LatticePropagatorSite> wl_eta_sliceSum, ws_eta_sliceSum;

          LatticePropagator wl_eta_Coulomb = env.toCoulombSink(gt, wl_eta);
          sliceSum(wl_eta_Coulomb, wl_eta_sliceSum, Tdir);
          wl_tEta_tEta = wl_eta_sliceSum[t_eta];

          LatticePropagator ws_eta_Coulomb = env.toCoulombSink(gt, ws_eta);
          sliceSum(ws_eta_Coulomb, ws_eta_sliceSum, Tdir);
          ws_tEta_tEta = ws_eta_sliceSum[t_eta];
        }

        for(const string &diagram: diagrams) {
          LatticeComplex rst_lat = diagram_funcs[diagram](Lxx, wl_K, ws_K, wl_eta, ws_eta, wl_tEta_tK, ws_tEta_tK, wl_tEta_tEta, ws_tEta_tEta);
          std::vector<LatticeComplexSite> rst_slice_sum;  // sum over time slice of xt
          sliceSum(rst_lat, rst_slice_sum, Tdir);
          for(int xt=0; xt<T; ++xt) table_4d[diagram][traj_idx][t_sep_idx][t_K][xt] = TensorRemove(rst_slice_sum[(t_K+xt)%T]); 
        }

      }
    }

    ++traj_idx;
  }

  // Note: In python code, add coefficents, sqrt(2)/2 or -sqrt(2)/2  !!!!! This code calculates only the trace, without coefficient
  // Note: table_4d[traj][t_sep][t_K][x_t]; already shifted -> x_t=0 is the position of kaon
  for(const string &diagram: diagrams) {
    std::cout << std::string(40, '*') << std::endl;
    std::cout << diagram << ": " << average_tK(table_4d[diagram]) << std::endl;
  }

  std::cout << GridLogMessage << "Finished!" << std::endl;

  Grid_finalize();


  return 0;
}
