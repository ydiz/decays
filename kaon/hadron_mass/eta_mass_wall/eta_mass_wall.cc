#include "../../kaon.h"

using namespace std;
using namespace Grid;



vector<int> get_available_trajs() {

  string path = "/global/cfs/cdirs/mp13/ydzhao/24ID/wall_l/";
	DIR *dir;
	dir = opendir(path.c_str());
	assert(dir!=NULL); // make sure directory exists
	struct dirent *entry;

	string subdir_name;
  vector<int> trajs;
	while ((entry = readdir (dir)) != NULL) {
		subdir_name = string(entry->d_name);
		if(std::isdigit(subdir_name[0])) trajs.push_back(std::stoi(subdir_name));
	}
	closedir(dir);

  sort(trajs.begin(), trajs.end());
  return trajs;
}

void load_int_cmd_argument(int argc, char **argv, const string &arg_name, int &var) {
  if( GridCmdOptionExists(argv, argv+argc, arg_name) ) {
    string arg = GridCmdOptionPayload(argv, argv+argc, arg_name);
    GridCmdOptionInt(arg, var);
  }
  else {
    std::cout << arg_name << " not specified; exiting" << std::endl;
    exit(0);
  }
}


int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

#ifdef CUTH_FREE_FIELD
  Env env("FreeField_8nt8");
#else
  Env env("24ID");
#endif

	// int traj_start = 2300, traj_end = 2300, traj_sep = 100; 
	// int traj_start = 2260, traj_end = 2640, traj_sep = 10; 
	int traj_start, traj_end, traj_sep = 10; 
  load_int_cmd_argument(argc, argv, "--traj_start", traj_start);
  load_int_cmd_argument(argc, argv, "--traj_end", traj_end);


	cout << string(20, '*') << endl;
	cout << "traj_start: " << traj_start << endl;
	cout << "traj_end: " << traj_end << endl;
	cout << "traj_sep: " << traj_sep << endl;
	cout << string(20, '*') << endl;

  const int T = env.lat_size[3];

  vector<int> available_trajs = get_available_trajs();
  std::cout << "Available trajectories: " << available_trajs << std::endl;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    if(find(available_trajs.begin(), available_trajs.end(), traj) == available_trajs.end()) continue;
       
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l', true);  // true: sink should be in Coulomb gauge
    std::vector<LatticePropagator> ws = env.get_wall('s', true);  // true: sink should be in Coulomb gauge

    vector<vector<LatticePropagatorSite>> wl_sliceSum(T), ws_sliceSum(T); // wl_sliceSum[t1][t2] is the propagator from t1 to t2, i.e. L(t2, t1)
    for(int t=0; t<T; ++t) {
      sliceSum(wl[t], wl_sliceSum[t], Tdir);
      sliceSum(ws[t], ws_sliceSum[t], Tdir);
    }

    vector<vector<double>> all_corr;   // all_corr[t_eta][t_sep], the second index is tsep, not position of the second pion
    for(int t_eta=0; t_eta<T; ++t_eta) { // iterate through the position of eta wall

      // vector<LatticePropagatorSite> wl_sliceSum, ws_sliceSum;
      // sliceSum(wl[t_eta], wl_sliceSum, Tdir);
      // sliceSum(ws[t_eta], ws_sliceSum, Tdir);

      vector<double> corr(T, 0.);
      for(int i=0; i<T; ++i) {  // position of another eta
        int sep = (i - t_eta + T) % T;
        LatticeComplexSite amp_connected = - 1./3. * trace(wl_sliceSum[t_eta][i] * adj(wl_sliceSum[t_eta][i])) - 2./3. * trace(ws_sliceSum[t_eta][i] * adj(ws_sliceSum[t_eta][i]));
        LatticeComplexSite amp_disconnected = 2./3. * ( trace(g5 * wl_sliceSum[t_eta][t_eta]) * trace(g5 * wl_sliceSum[i][i]) + trace(g5 * ws_sliceSum[t_eta][t_eta]) * trace(g5 * ws_sliceSum[i][i]) - trace(g5 * wl_sliceSum[t_eta][t_eta]) * trace(g5 * ws_sliceSum[i][i]) - trace(g5 * ws_sliceSum[t_eta][t_eta]) * trace(g5 * wl_sliceSum[i][i])  );
        corr[sep] = ( amp_connected()()() + amp_disconnected()()() ).real();
      }
      // std::cout << "t_eta = "<< t_eta << ": "  << corr << std::endl;
      all_corr.push_back(corr);
    }

    // average over all positions of Kaon
    vector<double> avg_corr(T, 0.);  // averaged over all positions of t_eta
    for(int t_sep=0; t_sep<T; ++t_sep) {
      for(int t_eta=0; t_eta<T; ++t_eta) avg_corr[t_sep] += all_corr[t_eta][t_sep];
      avg_corr[t_sep] /= double(T);
    }

    std::cout << GridLogMessage << "traj [" << traj << "]: " << avg_corr << std::endl;
  }

  std::cout << "Finished!" << std::endl;
  Grid_finalize();
  return 0;
}
