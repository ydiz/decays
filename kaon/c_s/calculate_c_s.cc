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
  Grid_init(&argc, &argv);

  // int traj_start = 2300, traj_end = 2400, traj_sep = 100; // for 24ID, kaon wall
  int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  int t_sep_min = 20, t_sep_max= 20;   // separation of pion and kaon walls
  std::cout << "t_sep_min: " << t_sep_min << std::endl;
  std::cout << "t_sep_max: " << t_sep_max << std::endl;
  std::cout << "t_sep_num: " << t_sep_max - t_sep_min + 1 << std::endl;

  Env env(gcoor, "24ID");
  // init_para(argc, argv, env);
  const int T = env.grid->_fdimensions[3];

  vector<vector<vector<vector<Complex>>>> diagram1_Q1, diagram1_Q2, diagram2_Q1, diagram2_Q2, diagram_sBar_d; // diagram1_Q1[traj][t_sep][t_K][x_t]; already shifted -> x_t=0 is the position of kaon
  resize_vec(diagram1_Q1, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram1_Q2, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram2_Q1, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram2_Q2, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram_sBar_d, traj_num, T, t_sep_min, t_sep_max); 

  int traj_idx = 0;
  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    // Read gauge transformation matrices for this trajectory
    LatticeColourMatrix gt = env.get_gaugeTransform();
    // LatticeColourMatrix gt(env.grid);
    // readScidac(gt, env.gauge_transform_path());

    LatticePropagator Lxx = env.get_Lxx();

    for(int t_K=0; t_K<T; ++t_K) { // iterate through the position of Kaon wall
      LatticePropagator wl_K = env.get_wall(t_K, 'l'); // L(x, t_K)
      LatticePropagator ws_K = env.get_wall(t_K, 's'); // H(x, t_K)

      // Calculate L(t_pi, t_K); the sink must be in Coulomb gauge
      LatticePropagator wl_K_Coulomb = env.toCoulombSink(gt, wl_K);
      std::vector<LatticePropagatorSite> wl_K_sliceSum;
      sliceSum(wl_K_Coulomb, wl_K_sliceSum, Tdir);

      for(int t_sep = t_sep_min; t_sep<=t_sep_max; ++t_sep) { // iterate through all possible separation between pion and kaon
        int t_pi = (t_K + t_sep) % T;

        LatticePropagator wl_pi = env.get_wall(t_pi, 'l'); // L(x, t_pi)
        LatticePropagatorSite wl_tpi_tK = wl_K_sliceSum[t_pi];  // L(t_pi, t_K)

        // diagram_sBar_d 
        LatticeComplex rst_sBar_d(env.grid);
        rst_sBar_d = trace(g5 * wl_tpi_tK * adj(ws_K) * g5 * wl_pi);

        std::vector<LatticeComplexSite> rst_sBar_d_slice_sum;
        sliceSum(rst_sBar_d, rst_sBar_d_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram_sBar_d[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_sBar_d_slice_sum[(t_K+xt)%T]); 

        // diagram1_Q1 and diagram1_Q2
        LatticePropagator tmp1 = wl_pi * adj(wl_pi);
        LatticePropagator tmp2 = wl_K * adj(ws_K);
        LatticeComplex rst_D1_Q1(env.grid), rst_D1_Q2(env.grid);
        rst_D1_Q1 = Zero(); rst_D1_Q2 = Zero();
        for(int mu=0; mu<4; ++mu) {
          rst_D1_Q1 += trace(gL[mu] * tmp1) * trace(gL[mu] * tmp2);
          rst_D1_Q2 += trace(gL[mu] * tmp1 * gL[mu] * tmp2);
        }

        std::vector<LatticeComplexSite> rst_D1_Q1_slice_sum, rst_D1_Q2_slice_sum;
        sliceSum(rst_D1_Q1, rst_D1_Q1_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram1_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D1_Q1_slice_sum[(t_K+xt)%T]); 

        sliceSum(rst_D1_Q2, rst_D1_Q2_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram1_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D1_Q2_slice_sum[(t_K+xt)%T]); 

        // diagram2_Q1 and diagram2_Q2

        LatticePropagator tmp3 = wl_pi * g5 * wl_tpi_tK * adj(ws_K);

        LatticeComplex rst_D2_Q1(env.grid), rst_D2_Q2(env.grid);
        rst_D2_Q1 = Zero(); rst_D2_Q2 = Zero();
        for(int mu=0; mu<4; ++mu) {
          rst_D2_Q1 += trace(gL[mu] * Lxx) * trace(gL[mu] * tmp3);
          rst_D2_Q2 += trace(gL[mu] * Lxx * gL[mu] * tmp3);
        }

        std::vector<LatticeComplexSite> rst_D2_Q1_slice_sum, rst_D2_Q2_slice_sum;
        sliceSum(rst_D2_Q1, rst_D2_Q1_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram2_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D2_Q1_slice_sum[(t_K+xt)%T]); 

        sliceSum(rst_D2_Q2, rst_D2_Q2_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram2_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_D2_Q2_slice_sum[(t_K+xt)%T]); 
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
  std::cout << "diagram sBar_d: " << diagram_sBar_d << std::endl;
  std::cout << std::string(40, '*') << std::endl;
  std::cout << "Finished!" << std::endl;

  Grid_finalize();

  return 0;
}
