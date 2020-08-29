#include "../kaon.h"

// using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;

std::vector<int> gcoor({24, 24, 24, 64});

void resize_vec(std::vector<std::vector<std::vector<std::vector<Grid::Complex>>>> &diagram, int traj_num, int T, int t_sep_min, int t_sep_max) {
  diagram.resize(traj_num);
  for(auto &x: diagram) {
    x.resize(t_sep_max - t_sep_min + 1);
    for(auto &y: x) {
      y.resize(T);
      for(auto &z: y) z.resize(T);
    }
  }
}


int main(int argc, char* argv[])
{
  // Grid_init(&argc, &argv);
  zyd_init_Grid_Qlattice(argc, argv);

  // int traj_start = 2300, traj_end = 2400, traj_sep = 100; // for 24ID, kaon wall
  int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  int t_sep_min = 20, t_sep_max= 20;   // separation of eta and kaon walls
  std::cout << "t_sep_min: " << t_sep_min << std::endl;
  std::cout << "t_sep_max: " << t_sep_max << std::endl;
  std::cout << "t_sep_num: " << t_sep_max - t_sep_min + 1 << std::endl;

  Env env("24ID");
  // init_para(argc, argv, env);
  const int T = env.grid->_fdimensions[3];

  vector<vector<vector<vector<Complex>>>> diagram1_Q1, diagram1_Q2, diagram2_Q1, diagram2_Q2, diagram3_Q1, diagram3_Q2, diagram4_Q1, diagram4_Q2, diagram_sBar_d_D1, diagram_sBar_d_D2, diagram_sBar_d_D3D4; // diagram1_Q1[traj][t_sep][t_K][x_t]; already shifted -> x_t=0 is the position of kaon
  resize_vec(diagram1_Q1, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram1_Q2, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram2_Q1, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram2_Q2, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram3_Q1, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram3_Q2, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram4_Q1, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram4_Q2, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram_sBar_d_D1, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram_sBar_d_D2, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram_sBar_d_D3D4, traj_num, T, t_sep_min, t_sep_max); 

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

      // Calculate L(t_eta, t_K); the sink must be in Coulomb gauge
      std::vector<LatticePropagatorSite> ws_K_sliceSum;
      {
        LatticePropagator ws_K_Coulomb = env.toCoulombSink(gt, ws_K);
        sliceSum(ws_K_Coulomb, ws_K_sliceSum, Tdir);
      }

      for(int t_sep = t_sep_min; t_sep<=t_sep_max; ++t_sep) { // iterate through all possible separation between eta and kaon
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


        // diagram_sBar_d_D1 and diagram_sBar_d_D2
        LatticeComplex rst_sBar_d_D1(env.grid), rst_sBar_d_D2(env.grid);
        rst_sBar_d_D1 = trace(g5 * wl_tEta_tK * adj(ws_K) * g5 * wl_eta);
        rst_sBar_d_D2 = trace(g5 * adj(ws_eta) * g5 * wl_K * adj(ws_tEta_tK));

        std::vector<LatticeComplexSite> rst_sBar_d_D1_slice_sum, rst_sBar_d_D2_slice_sum;
        sliceSum(rst_sBar_d_D1, rst_sBar_d_D1_slice_sum, Tdir);
        sliceSum(rst_sBar_d_D2, rst_sBar_d_D2_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram_sBar_d_D1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_sBar_d_D1_slice_sum[(t_K+xt)%T]); 
        for(int xt=0; xt<T; ++xt) diagram_sBar_d_D2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_sBar_d_D2_slice_sum[(t_K+xt)%T]); 

        // diagram_sBar_d_D3D4
        LatticeComplex rst_sBar_d_D3D4(env.grid);
        rst_sBar_d_D3D4 =  trace(g5 * (wl_tEta_tEta - ws_tEta_tEta) ) * trace(g5 * wl_K * adj(ws_K));

        std::vector<LatticeComplexSite> rst_sBar_d_D3D4_slice_sum;
        sliceSum(rst_sBar_d_D3D4, rst_sBar_d_D3D4_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram_sBar_d_D3D4[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_sBar_d_D3D4_slice_sum[(t_K+xt)%T]); 


        // diagram1_Q1 and diagram1_Q2
        LatticePropagator tmp1 = wl_eta * adj(wl_eta);
        LatticePropagator tmp2 = wl_K * adj(ws_K);
        LatticeComplex rst_D1_Q1(env.grid), rst_D1_Q2(env.grid);
        rst_D1_Q1 = Zero(); rst_D1_Q2 = Zero();
        for(int mu=0; mu<4; ++mu) {
          rst_D1_Q1 += trace(gL[mu] * tmp1) * trace(gL[mu] * tmp2);
          rst_D1_Q2 += trace(gL[mu] * tmp1 * gL[mu] * tmp2);
        }

        std::vector<LatticeComplexSite> rst_D1_Q1_slice_sum, rst_D1_Q2_slice_sum;
        sliceSum(rst_D1_Q1, rst_D1_Q1_slice_sum, Tdir);
        sliceSum(rst_D1_Q2, rst_D1_Q2_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram1_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D1_Q1_slice_sum[(t_K+xt)%T]); 
        for(int xt=0; xt<T; ++xt) diagram1_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D1_Q2_slice_sum[(t_K+xt)%T]); 

        // diagram2_Q1 and diagram2_Q2

        LatticePropagator tmp3 = wl_eta * g5 * wl_tEta_tK * adj(ws_K);

        LatticeComplex rst_D2_Q1(env.grid), rst_D2_Q2(env.grid);
        rst_D2_Q1 = Zero(); rst_D2_Q2 = Zero();
        for(int mu=0; mu<4; ++mu) {
          rst_D2_Q1 += trace(gL[mu] * Lxx) * trace(gL[mu] * tmp3);
          rst_D2_Q2 += trace(gL[mu] * Lxx * gL[mu] * tmp3);
        }

        std::vector<LatticeComplexSite> rst_D2_Q1_slice_sum, rst_D2_Q2_slice_sum;
        sliceSum(rst_D2_Q1, rst_D2_Q1_slice_sum, Tdir);
        sliceSum(rst_D2_Q2, rst_D2_Q2_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram2_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D2_Q1_slice_sum[(t_K+xt)%T]); 
        for(int xt=0; xt<T; ++xt) diagram2_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D2_Q2_slice_sum[(t_K+xt)%T]); 

        // diagram3_Q1 and diagram3_Q2

        LatticePropagator tmp4 = wl_K * adj(ws_tEta_tK) * g5 * adj(ws_eta);

        LatticeComplex rst_D3_Q1(env.grid), rst_D3_Q2(env.grid);
        rst_D3_Q1 = Zero(); rst_D3_Q2 = Zero();
        for(int mu=0; mu<4; ++mu) {
          rst_D3_Q1 += trace(gL[mu] * Lxx) * trace(gL[mu] * tmp4);
          rst_D3_Q2 += trace(gL[mu] * Lxx * gL[mu] * tmp4);
        }

        std::vector<LatticeComplexSite> rst_D3_Q1_slice_sum, rst_D3_Q2_slice_sum;
        sliceSum(rst_D3_Q1, rst_D3_Q1_slice_sum, Tdir);
        sliceSum(rst_D3_Q2, rst_D3_Q2_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram3_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D3_Q1_slice_sum[(t_K+xt)%T]); 
        for(int xt=0; xt<T; ++xt) diagram3_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D3_Q2_slice_sum[(t_K+xt)%T]); 

        // diagram4_Q1 and diagram4_Q2

        LatticeComplex rst_D4_Q1(env.grid), rst_D4_Q2(env.grid);
        rst_D4_Q1 = Zero(); rst_D4_Q2 = Zero();
        for(int mu=0; mu<4; ++mu) {
          rst_D4_Q1 += trace(g5 * (wl_tEta_tEta - ws_tEta_tEta) ) * trace(gL[mu] * Lxx) * trace(gL[mu] * wl_K * adj(ws_K) );
          rst_D4_Q2 += trace(g5 * (wl_tEta_tEta - ws_tEta_tEta) ) * trace(gL[mu] * Lxx * gL[mu] * wl_K * adj(ws_K) );
        }

        std::vector<LatticeComplexSite> rst_D4_Q1_slice_sum, rst_D4_Q2_slice_sum;
        sliceSum(rst_D4_Q1, rst_D4_Q1_slice_sum, Tdir);
        sliceSum(rst_D4_Q2, rst_D4_Q2_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram4_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D4_Q1_slice_sum[(t_K+xt)%T]); 
        for(int xt=0; xt<T; ++xt) diagram4_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D4_Q2_slice_sum[(t_K+xt)%T]); 
      }
    }

    ++traj_idx;
  }

// Note: In python code, add coefficents, sqrt(2)/2 or -sqrt(2)/2  !!!!! This code calculates only the trace, without coefficient
// diagram1_Q1[traj][t_sep][t_K][x_t]; already shifted -> x_t=0 is the position of kaon
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram1 Q1: " << diagram1_Q1 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram1 Q2: " << diagram1_Q2 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram2 Q1: " << diagram2_Q1 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram2 Q2: " << diagram2_Q2 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram3 Q1: " << diagram3_Q1 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram3 Q2: " << diagram3_Q2 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram4 Q1: " << diagram4_Q1 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram4 Q2: " << diagram4_Q2 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram sBar_d_D1: " << diagram_sBar_d_D1 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram sBar_d_D2: " << diagram_sBar_d_D2 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "diagram sBar_d_D3D4: " << diagram_sBar_d_D3D4 << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "Finished!" << std::endl;

  Grid_finalize();

  return 0;
}
