// #include <headers/headers.h>
//
// namespace Grid {
// namespace QCD {
//
// std::vector<double> kaon_corr(int traj, GridBase *grid) {
//   int T = grid->_fdimensions[Tdir];
// 	LatticePropagator prop_d(grid);
// 	LatticePropagator prop_s(grid);
//
//   std::vector<double> corr(T, 0.);
//
//   for(int t=0; t<T; ++t) {
//
//     // std::string path_d = wall_path_ud_32D(traj, t);
//     // std::string path_s = wall_path_strange_32D(traj, t);
//     std::string path_d = wall_path_l_24ID(traj, t);
//     std::string path_s = wall_path_s_24ID(traj, t);
//     read_qlat_propagator(prop_d, path_d);
//     read_qlat_propagator(prop_s, path_s);
//
//     std::vector<typename LatticePropagator::vector_object::scalar_object> slice_sum_d;
//     sliceSum(prop_d, slice_sum_d, Tdir);
//     std::vector<typename LatticePropagator::vector_object::scalar_object> slice_sum_s;
//     sliceSum(prop_s, slice_sum_s, Tdir);
//
//     for(int i=0; i<T; ++i) {
//       int sep = (i < t) ? i - t + T : i - t;
//       double tmp = TensorRemove(trace(slice_sum_d[i] * adj(slice_sum_s[i]))).real();
//       corr[sep] += tmp;
//     }
//   }
//
//   for(auto &x: corr) x /= T;
//   std::cout << "Wall to Wall Correlator[traj=" << traj << "]:"<< std::endl;
//   std::cout << corr << std::endl;
//
//   return corr;
// }
//
//
// std::vector<double> pion_corr(int traj, GridBase *grid) {
//   int T = grid->_fdimensions[Tdir];
// 	LatticePropagator prop(grid);
//
//   std::vector<double> corr(T, 0.);
//
//   for(int t=0; t<T; ++t) {
//
//     // std::string path = wall_path_ud_32D(traj, t);
//     // std::string path = wall_path_24ID(traj, t);
//     std::string path = wall_path_ud_48I(traj, t);
//     read_qlat_propagator(prop, path);
//
//     std::vector<typename LatticePropagator::vector_object::scalar_object> slice_sum;
//     sliceSum(prop, slice_sum, Tdir);
//
//     for(int i=0; i<T; ++i) { // slice_sum[i] = P(i, t)
//       int sep = (i < t) ? i - t + T : i - t;
//       double tmp = TensorRemove(trace(slice_sum[i] * adj(slice_sum[i]))).real();
//       corr[sep] += tmp;
//     }
//   }
//
//   for(auto &x: corr) x /= T;
//   std::cout << "Wall to Wall Correlator[traj=" << traj << "]:"<< std::endl;
//   std::cout << corr << std::endl;
//
//   return corr;
// }
//
//
// }}


#include "kaon.h"

// using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;



int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  vector<int> gcoor({24, 24, 24, 64});
	int traj_start = 2300, traj_end = 2300, traj_sep = 100; 
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

	cout << string(20, '*') << endl;
	cout << "traj_start: " << traj_start << endl;
	cout << "traj_end: " << traj_end << endl;
	cout << "traj_sep: " << traj_sep << endl;
	cout << "traj_num: " << traj_num << endl;
	cout << string(20, '*') << endl;

  Env env(gcoor, "24ID");
  // init_para(argc, argv, env);
  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    vector<vector<double>> all_corr;   // all_corr[t_K][t_sep]
    env.setup_traj(traj);

    LatticeColourMatrix gt(env.grid);
    readScidac(gt, env.gauge_transform_path());

    for(int t_K=0; t_K<T; ++t_K) { // iterate through the position of Kaon wall
      LatticePropagator wl_K = env.get_wall(t_K, 'l'); // L(x, t_K)
      LatticePropagator ws_K = env.get_wall(t_K, 's'); // H(x, t_K)

      vector<LatticePropagatorSite> wl_K_sliceSum, ws_K_sliceSum;
      LatticePropagator wl_K_Coulomb = env.toCoulombSink(gt, wl_K);
      sliceSum(wl_K_Coulomb, wl_K_sliceSum, Tdir);
      LatticePropagator ws_K_Coulomb = env.toCoulombSink(gt, ws_K);
      sliceSum(ws_K_Coulomb, ws_K_sliceSum, Tdir);

      vector<double> corr(T, 0.);
      for(int i=0; i<T; ++i) {
        int sep = (i < t_K) ? i - t_K + T : i - t_K;
        double tmp = TensorRemove(trace(ws_K_sliceSum[i] * adj(wl_K_sliceSum[i]))).real();
        corr[sep] = tmp;
      }
      all_corr.push_back(corr);

    }

    // average over all positions of Kaon
    vector<double> avg_corr(T, 0.);  // averaged over all positions of t_K
    for(int t_sep=0; t_sep<T; ++t_sep) {
      for(int t_K=0; t_K<T; ++t_K) avg_corr[t_sep] += all_corr[t_K][t_sep];
      avg_corr[t_sep] /= double(T);
    }

    std::cout << GridLogMessage << "traj [" << traj << "]: " << avg_corr << std::endl;
  }

	// cout << string(30, '*') << endl;
	// cout << "wall to wall correlators over "<< traj_num << " trajectoies:" << endl;
	// cout << all_corr << endl;


  std::cout << "Finished!" << std::endl;
  Grid_finalize();
  return 0;
}
