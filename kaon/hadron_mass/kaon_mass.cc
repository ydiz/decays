// On 16 nodes, needs ~8 hours to run


#include "../kaon.h"

using namespace std;
using namespace Grid;


int main(int argc, char* argv[])
{
  // zyd_init_Grid_Qlattice(argc, argv);
  Grid_init(&argc, &argv);

  vector<int> trajs = {1010, 1030, 1040, 1050, 1070, 1080, 1090, 1110, 1120, 1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1260, 1270, 1280, 1290, 1300, 1310, 1320, 1330, 1350, 1360, 1370, 1380, 1390, 1400};
  // vector<int> trajs = {1010, 1030, 1040, 1050, 1070, 1080, 1090, 1110, 1120, 1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1260, 1270, 1280, 1290, 1300, 1310, 1320, 1330, 1350, 1360, 1370, 1380, 1390, 1400, 1410, 1420, 1430, 1440, 1450, 1460, 1470, 1480, 1490, 1500, 1510, 1520, 1530, 1540, 1550, 1560, 1570, 1580, 1590, 1600, 1610, 1620, 1630, 1640, 1650, 1660, 1670, 1680, 1690, 1700, 1710, 1720, 1730, 1740, 1750, 1760, 1770, 1780, 1790, 1800, 1810, 1820, 1830, 1840, 1850, 1860, 1870, 1880, 1890, 1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100, 2110, 2120, 2130, 2140, 2150, 2160, 2170, 2180, 2190, 2200, 2210, 2220, 2230, 2240, 2250};

  std::cout << "trajs: " << trajs << std::endl;
	cout << "traj_num: " << trajs.size() << endl;

  Env env("24ID");
  const int T = env.lat_size[3];

  // for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
  for(int traj: trajs) {
    // if(find(traj_skips.begin(), traj_skips.end(), traj) != traj_skips.end()) continue;

    vector<vector<double>> all_corr;   // all_corr[t_K][t_sep], the second index is tsep, not position of the second kaon
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l', true);  // true: sink should be in Coulomb gauge
    std::vector<LatticePropagator> ws = env.get_wall('s', true);

    for(int t_K=0; t_K<T; ++t_K) { // iterate through the position of Kaon wall

      vector<LatticePropagatorSite> wl_K_sliceSum, ws_K_sliceSum;
      sliceSum(wl[t_K], wl_K_sliceSum, Tdir);
      sliceSum(ws[t_K], ws_K_sliceSum, Tdir);

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

  std::cout << "Finished!" << std::endl;
  Grid_finalize();
  return 0;
}
