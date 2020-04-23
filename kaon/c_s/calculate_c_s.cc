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

// std::vector<std::vector<std::vector<Grid::Complex>>> averaging_over_tK(std::vector<std::vector<std::vector<std::vector<Grid::Complex>>>> &diagram) {
//   assert(diagram.size() !=0);
//   int T = diagram[0][0][0].size();
//
//   std::vector<std::vector<std::vector<Grid::Complex>>> rst;
//   rst.resize(diagram.size());
//   for(i=0; i<rst.size(); ++i) {
//     rst[i].resize(diagram[i].size());
//
//     for(int sep_idx=0; sep_idx<rst[i].size(); ++sep_idx) {
//
//       rst[i][sep_idx].resize(T, 0.);
//
//       for(int tK=0; tK<T; ++tK) 
//         for(int xt=0; xt<T; ++xt)
//           rst[i][sep_idx][xt] += diagram[i][sep_idx][tK][xt];
//
//       for(int xt=0; xt<T; ++xt) rst[i][sep_idx][xt] /= double(T);
//     }
//   }
//   return rst;
// }


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

  vector<vector<vector<vector<Complex>>>> diagram1_Q1, diagram1_Q2, diagram_sBar_d; // diagram1_Q1[traj][t_sep][t_K][x_t]; already shifted -> x_t=0 is the position of kaon
  resize_vec(diagram1_Q1, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram1_Q2, traj_num, T, t_sep_min, t_sep_max); 
  resize_vec(diagram_sBar_d, traj_num, T, t_sep_min, t_sep_max); 

  int traj_idx = 0;
  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    // std::vector<LatticePropagator> wl = env.get_wall('l');
    // std::cout << "xxxx" << std::endl;
    // std::vector<LatticePropagator> ws = env.get_wall('s'); // will have error when creating the second vector<LatticePropagator>; Do not know why.
    // std::cout << "xxxx" << std::endl;
    // assert(0);

    LatticeColourMatrix gt(env.grid);
    readScidac(gt, env.gauge_transform_path());

    for(int t_K=0; t_K<T; ++t_K) { // iterate through the position of Kaon wall
    // for(int t_K=0; t_K<1; ++t_K) {
      LatticePropagator wl_K = env.get_wall(t_K, 'l');
      LatticePropagator ws_K = env.get_wall(t_K, 's');

      for(int t_sep = t_sep_min; t_sep<=t_sep_max; ++t_sep) { // iterate through all possible separation between pion and kaon
        int t_pi = (t_K + t_sep) % T;

        LatticePropagator wl_pi = env.get_wall(t_pi, 'l');

        // diagram1_Q1 and diagram1_Q2
        LatticePropagator tmp1 = wl_pi * adj(wl_pi);
        LatticePropagator tmp2 = wl_K * adj(ws_K);
        LatticeComplex rst_Q1(env.grid), rst_Q2(env.grid);
        rst_Q1 = Zero(); rst_Q2 = Zero();
        for(int mu=0; mu<4; ++mu) {
          rst_Q1 += trace(gL[mu] * tmp1) * trace(gL[mu] * tmp2);
          rst_Q2 += trace(gL[mu] * tmp1 * gL[mu] * tmp2);
        }
        std::vector<typename LatticeComplex::vector_object::scalar_object> rst_Q1_slice_sum, rst_Q2_slice_sum;
        sliceSum(rst_Q1, rst_Q1_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram1_Q1[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_Q1_slice_sum[(t_K+xt)%T]); 
        sliceSum(rst_Q2, rst_Q2_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram1_Q2[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_Q2_slice_sum[(t_K+xt)%T]); 

        // diagram_sBar_d 

        // For L(t_pi, t_K), the 
        LatticePropagator wl_K_Coulomb = env.toCoulombSink(gt, wl_K);

        std::vector<typename LatticePropagator::vector_object::scalar_object> wl_K_sliceSum;
        // sliceSum(wl_K, wl_K_sliceSum, Tdir);
        sliceSum(wl_K_Coulomb, wl_K_sliceSum, Tdir);
        typename LatticePropagator::vector_object::scalar_object wl_tpi_tK = wl_K_sliceSum[t_pi];
        
        LatticeComplex rst_sBar_d(env.grid);
        rst_sBar_d = trace(g5 * wl_tpi_tK * adj(ws_K) * g5 * wl_pi);
        std::vector<typename LatticeComplex::vector_object::scalar_object> rst_sBar_d_slice_sum;
        sliceSum(rst_sBar_d, rst_sBar_d_slice_sum, Tdir);
        for(int xt=0; xt<T; ++xt) diagram_sBar_d[traj_idx][t_sep - t_sep_min][t_K][xt] = TensorRemove(rst_sBar_d_slice_sum[(t_K+xt)%T]); 
      }
    }

    ++traj_idx;
  }
    // It might be better to add coefficients in python code, not here. // so that in here we calculate only contractions.
// TBD: add coefficents, sqrt(2) / 2 or - sqrt(2) / 2  !!!!!
// diagram1_Q1[traj][t_sep][t_K][x_t]; already shifted -> x_t=0 is the position of kaon
  std::cout << "diagram1 Q1: " << diagram1_Q1 << std::endl;
  std::cout << "diagram1 Q2: " << diagram1_Q2 << std::endl;
  std::cout << "diagram sBar_d: " << diagram_sBar_d << std::endl;
  std::cout << "Finished!" << std::endl;

  Grid_finalize();

  // std::this_thread::sleep_for(std::chrono::seconds(1));  // sleep for 1 second

  return 0;
}
