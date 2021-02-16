#include "../kaon.h"

using namespace std;
using namespace Grid;

int main(int argc, char* argv[])
{
  zyd_init_Grid_Qlattice(argc, argv);

	// int traj_start = 2300, traj_end = 2300, traj_sep = 100; 
	int traj_start = 2260, traj_end = 2640, traj_sep = 10; 
  std::vector<int> traj_skips = {2360, 2420, 2520, 2540, 2580};   // 2420 has no strange wall source
  int traj_num = (traj_end - traj_start) / traj_sep + 1 - traj_skips.size();

	cout << string(20, '*') << endl;
	cout << "traj_start: " << traj_start << endl;
	cout << "traj_end: " << traj_end << endl;
	cout << "traj_sep: " << traj_sep << endl;
	cout << "traj_skip: " << traj_skips << endl;
	cout << "traj_num: " << traj_num << endl;
	cout << string(20, '*') << endl;

  Env env("24ID");
  const int T = env.lat_size[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    if(find(traj_skips.begin(), traj_skips.end(), traj) != traj_skips.end()) continue;
       
    vector<vector<double>> all_corr;   // all_corr[t_K][t_sep], the second index is tsep, not position of the second pion
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l', true);  // true: sink should be in Coulomb gauge

    for(int t_K=0; t_K<T; ++t_K) { // iterate through the position of Kaon wall

      vector<LatticePropagatorSite> wl_sliceSum;
      sliceSum(wl[t_K], wl_sliceSum, Tdir);

      vector<double> corr(T, 0.);
      for(int i=0; i<T; ++i) {
        int sep = (i < t_K) ? i - t_K + T : i - t_K;
        double tmp = TensorRemove(trace(wl_sliceSum[i] * adj(wl_sliceSum[i]))).real();
        corr[sep] = tmp;
      }
      // std::cout << "t_K = "<< t_K << ": "  << corr << std::endl;
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

  std::cout << "Finished!" << std::endl;
  Grid_finalize();
  return 0;
}
