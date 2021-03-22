#include "../kaon.h"
#include "c_s_eta_diagrams.h"

using namespace std;
using namespace Grid;

// void resize_vec(std::vector<std::vector<std::vector<std::vector<Grid::Complex>>>> &diagram, int traj_num, int T, int t_sep_min, int t_sep_max) {
//   diagram.resize(traj_num);
//   for(auto &x: diagram) {
//     x.resize(t_sep_max - t_sep_min + 1);
//     for(auto &y: x) {
//       y.resize(T);
//       for(auto &z: y) z.resize(T);
//     }
//   }
// }

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
  // int t_sep_min = 1, t_sep_max= 3; 
  vector<int> t_seps = {1,2,3}; 
  Env env("FreeField_8nt8");
#else
  // int t_sep_min = 10, t_sep_max= 24; // For eta: 10-18; for pion 18-24  // separation of eta and kaon walls
  vector<int> t_seps = {10, 12, 14, 16, 18, 20, 22, 24};  // For eta: 10-18; for pion 18-24  // separation of eta and kaon walls
  Env env("24ID");
#endif

  const int T = env.grid->_fdimensions[3];

  int N_t_seps = t_seps.size();
  std::cout << "t_seps: " << t_seps << std::endl;
  // std::cout << "t_sep_min: " << t_sep_min << std::endl;
  // std::cout << "t_sep_max: " << t_sep_max << std::endl;
  // std::cout << "t_sep_num: " << t_sep_max - t_sep_min + 1 << std::endl;

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


  // vector<vector<vector<vector<Complex>>>> diagram1_Q1, diagram1_Q2, diagram2_Q1, diagram2_Q2, diagram3_Q1, diagram3_Q2, diagram4_Q1, diagram4_Q2, diagram_sBar_d_D1, diagram_sBar_d_D2, diagram_sBar_d_D3D4; // diagram1_Q1[traj][t_sep][t_K][x_t]; already shifted -> x_t=0 is the position of kaon
  // resize_vec(diagram1_Q1, traj_num, T, t_sep_min, t_sep_max); 
  // resize_vec(diagram1_Q2, traj_num, T, t_sep_min, t_sep_max); 
  // resize_vec(diagram2_Q1, traj_num, T, t_sep_min, t_sep_max); 
  // resize_vec(diagram2_Q2, traj_num, T, t_sep_min, t_sep_max); 
  // resize_vec(diagram3_Q1, traj_num, T, t_sep_min, t_sep_max); 
  // resize_vec(diagram3_Q2, traj_num, T, t_sep_min, t_sep_max); 
  // resize_vec(diagram4_Q1, traj_num, T, t_sep_min, t_sep_max); 
  // resize_vec(diagram4_Q2, traj_num, T, t_sep_min, t_sep_max); 
  // resize_vec(diagram_sBar_d_D1, traj_num, T, t_sep_min, t_sep_max); 
  // resize_vec(diagram_sBar_d_D2, traj_num, T, t_sep_min, t_sep_max); 
  // resize_vec(diagram_sBar_d_D3D4, traj_num, T, t_sep_min, t_sep_max); 

  int traj_idx = 0;
  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    // Read gauge transformation matrices for this trajectory
    LatticeColourMatrix gt = env.get_gaugeTransform();
    LatticePropagator Lxx = env.get_Lxx();

    vector<LatticePropagator> wl = env.get_wall('l');
    vector<LatticePropagator> ws = env.get_wall('s');


    for(int t_K=0; t_K<T; ++t_K) { // iterate through the position of Kaon wall
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

        // // diagram_sBar_d_D1 and diagram_sBar_d_D2
        // LatticeComplex rst_sBar_d_D1(env.grid), rst_sBar_d_D2(env.grid);
        // rst_sBar_d_D1 = trace(g5 * wl_tEta_tK * adj(ws_K) * g5 * wl_eta);
        // rst_sBar_d_D2 = trace(g5 * adj(ws_eta) * g5 * wl_K * adj(ws_tEta_tK));
        //
        // std::vector<LatticeComplexSite> rst_sBar_d_D1_slice_sum, rst_sBar_d_D2_slice_sum;
        // sliceSum(rst_sBar_d_D1, rst_sBar_d_D1_slice_sum, Tdir);
        // sliceSum(rst_sBar_d_D2, rst_sBar_d_D2_slice_sum, Tdir);
        // for(int xt=0; xt<T; ++xt) diagram_sBar_d_D1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_sBar_d_D1_slice_sum[(t_K+xt)%T]); 
        // for(int xt=0; xt<T; ++xt) diagram_sBar_d_D2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_sBar_d_D2_slice_sum[(t_K+xt)%T]); 
        //
        // // diagram_sBar_d_D3D4
        // LatticeComplex rst_sBar_d_D3D4(env.grid);
        // rst_sBar_d_D3D4 =  trace(g5 * (wl_tEta_tEta - ws_tEta_tEta) ) * trace(g5 * wl_K * adj(ws_K));
        //
        // std::vector<LatticeComplexSite> rst_sBar_d_D3D4_slice_sum;
        // sliceSum(rst_sBar_d_D3D4, rst_sBar_d_D3D4_slice_sum, Tdir);
        // for(int xt=0; xt<T; ++xt) diagram_sBar_d_D3D4[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_sBar_d_D3D4_slice_sum[(t_K+xt)%T]); 
        //
        //
        // // diagram1_Q1 and diagram1_Q2
        // LatticePropagator tmp1 = wl_eta * adj(wl_eta);
        // LatticePropagator tmp2 = wl_K * adj(ws_K);
        // LatticeComplex rst_D1_Q1(env.grid), rst_D1_Q2(env.grid);
        // rst_D1_Q1 = Zero(); rst_D1_Q2 = Zero();
        // for(int mu=0; mu<4; ++mu) {
        //   rst_D1_Q1 += trace(gL[mu] * tmp1) * trace(gL[mu] * tmp2);
        //   rst_D1_Q2 += trace(gL[mu] * tmp1 * gL[mu] * tmp2);
        // }
        //
        // std::vector<LatticeComplexSite> rst_D1_Q1_slice_sum, rst_D1_Q2_slice_sum;
        // sliceSum(rst_D1_Q1, rst_D1_Q1_slice_sum, Tdir);
        // sliceSum(rst_D1_Q2, rst_D1_Q2_slice_sum, Tdir);
        // for(int xt=0; xt<T; ++xt) diagram1_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D1_Q1_slice_sum[(t_K+xt)%T]); 
        // for(int xt=0; xt<T; ++xt) diagram1_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D1_Q2_slice_sum[(t_K+xt)%T]); 
        //
        // // diagram2_Q1 and diagram2_Q2
        //
        // LatticePropagator tmp3 = wl_eta * g5 * wl_tEta_tK * adj(ws_K);
        //
        // LatticeComplex rst_D2_Q1(env.grid), rst_D2_Q2(env.grid);
        // rst_D2_Q1 = Zero(); rst_D2_Q2 = Zero();
        // for(int mu=0; mu<4; ++mu) {
        //   rst_D2_Q1 += trace(gL[mu] * Lxx) * trace(gL[mu] * tmp3);
        //   rst_D2_Q2 += trace(gL[mu] * Lxx * gL[mu] * tmp3);
        // }
        //
        // std::vector<LatticeComplexSite> rst_D2_Q1_slice_sum, rst_D2_Q2_slice_sum;
        // sliceSum(rst_D2_Q1, rst_D2_Q1_slice_sum, Tdir);
        // sliceSum(rst_D2_Q2, rst_D2_Q2_slice_sum, Tdir);
        // for(int xt=0; xt<T; ++xt) diagram2_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D2_Q1_slice_sum[(t_K+xt)%T]); 
        // for(int xt=0; xt<T; ++xt) diagram2_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D2_Q2_slice_sum[(t_K+xt)%T]); 
        //
        // // diagram3_Q1 and diagram3_Q2
        //
        // LatticePropagator tmp4 = wl_K * adj(ws_tEta_tK) * g5 * adj(ws_eta);
        //
        // LatticeComplex rst_D3_Q1(env.grid), rst_D3_Q2(env.grid);
        // rst_D3_Q1 = Zero(); rst_D3_Q2 = Zero();
        // for(int mu=0; mu<4; ++mu) {
        //   rst_D3_Q1 += trace(gL[mu] * Lxx) * trace(gL[mu] * tmp4);
        //   rst_D3_Q2 += trace(gL[mu] * Lxx * gL[mu] * tmp4);
        // }
        //
        // std::vector<LatticeComplexSite> rst_D3_Q1_slice_sum, rst_D3_Q2_slice_sum;
        // sliceSum(rst_D3_Q1, rst_D3_Q1_slice_sum, Tdir);
        // sliceSum(rst_D3_Q2, rst_D3_Q2_slice_sum, Tdir);
        // for(int xt=0; xt<T; ++xt) diagram3_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D3_Q1_slice_sum[(t_K+xt)%T]); 
        // for(int xt=0; xt<T; ++xt) diagram3_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D3_Q2_slice_sum[(t_K+xt)%T]); 
        //
        // // diagram4_Q1 and diagram4_Q2
        //
        // LatticeComplex rst_D4_Q1(env.grid), rst_D4_Q2(env.grid);
        // rst_D4_Q1 = Zero(); rst_D4_Q2 = Zero();
        // for(int mu=0; mu<4; ++mu) {
        //   rst_D4_Q1 += trace(g5 * (wl_tEta_tEta - ws_tEta_tEta) ) * trace(gL[mu] * Lxx) * trace(gL[mu] * wl_K * adj(ws_K) );
        //   rst_D4_Q2 += trace(g5 * (wl_tEta_tEta - ws_tEta_tEta) ) * trace(gL[mu] * Lxx * gL[mu] * wl_K * adj(ws_K) );
        // }
        //
        // std::vector<LatticeComplexSite> rst_D4_Q1_slice_sum, rst_D4_Q2_slice_sum;
        // sliceSum(rst_D4_Q1, rst_D4_Q1_slice_sum, Tdir);
        // sliceSum(rst_D4_Q2, rst_D4_Q2_slice_sum, Tdir);
        // for(int xt=0; xt<T; ++xt) diagram4_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D4_Q1_slice_sum[(t_K+xt)%T]); 
        // for(int xt=0; xt<T; ++xt) diagram4_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D4_Q2_slice_sum[(t_K+xt)%T]); 
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
  // // std::cout << "diagram1 Q1: " << diagram1_Q1 << std::endl;
  // std::cout << std::string(40, '*') << std::endl;
  // std::cout << "diagram1 Q2: " << average_tK(diagram1_Q2) << std::endl;
  // std::cout << std::string(40, '*') << std::endl;
  // std::cout << "diagram2 Q1: " << average_tK(diagram2_Q1) << std::endl;
  // std::cout << std::string(40, '*') << std::endl;
  // std::cout << "diagram2 Q2: " << average_tK(diagram2_Q2) << std::endl;
  // std::cout << std::string(40, '*') << std::endl;
  // std::cout << "diagram3 Q1: " << average_tK(diagram3_Q1) << std::endl;
  // std::cout << std::string(40, '*') << std::endl;
  // std::cout << "diagram3 Q2: " << average_tK(diagram3_Q2) << std::endl;
  // std::cout << std::string(40, '*') << std::endl;
  // std::cout << "diagram4 Q1: " << average_tK(diagram4_Q1) << std::endl;
  // std::cout << std::string(40, '*') << std::endl;
  // std::cout << "diagram4 Q2: " << average_tK(diagram4_Q2) << std::endl;
  // std::cout << std::string(40, '*') << std::endl;
  // std::cout << "diagram sBar_d_D1: " << average_tK(diagram_sBar_d_D1) << std::endl;
  // std::cout << std::string(40, '*') << std::endl;
  // std::cout << "diagram sBar_d_D2: " << average_tK(diagram_sBar_d_D2) << std::endl;
  // std::cout << std::string(40, '*') << std::endl;
  // std::cout << "diagram sBar_d_D3D4: " << average_tK(diagram_sBar_d_D3D4)<< std::endl;
  // std::cout << std::string(40, '*') << std::endl;

  std::cout << GridLogMessage << "Finished!" << std::endl;

  Grid_finalize();


  return 0;
}
