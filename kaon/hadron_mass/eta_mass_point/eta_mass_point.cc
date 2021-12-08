#include "../../kaon.h"

using namespace std;
using namespace Grid;



// vector<int> get_available_trajs() {
//
//   string path = "/global/cfs/cdirs/mp13/ydzhao/24ID/point_l/";
// 	DIR *dir;
// 	dir = opendir(path.c_str());
// 	assert(dir!=NULL); // make sure directory exists
// 	struct dirent *entry;
//
// 	string subdir_name;
//   vector<int> trajs;
// 	while ((entry = readdir (dir)) != NULL) {
// 		subdir_name = string(entry->d_name);
// 		if(std::isdigit(subdir_name[0])) trajs.push_back(std::stoi(subdir_name));
// 	}
// 	closedir(dir);
//
//   sort(trajs.begin(), trajs.end());
//   return trajs;
// }

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

  assert(env.grid->_Nprocessors == 1); // Only allow one process to run; no MPI

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

  // vector<int> available_trajs = get_available_trajs();
  // std::cout << "Available trajectories: " << available_trajs << std::endl;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    // if(find(available_trajs.begin(), available_trajs.end(), traj) == available_trajs.end()) continue;
       
    env.setup_traj(traj);

    if(env.xgs_l.size() <= 500) {
      std::cout << "trajectory " << traj << " has only " << env.xgs_l.size() << " point sources; Skipping this trajectory" << std::endl;
      continue;
    }

    map<vector<int>, map<vector<int>, LatticePropagatorSite>> pl_site, ps_site; // pl_site[x][y] is L(y, x), the point source prop from x to y
    int num_pt_src = 0;
    for(const auto &x: env.xgs_l) {  
      std::cout << "# Point source: " << num_pt_src  << " Point source: " << x << std::endl;
      LatticePropagator pl = env.get_point(x, 'l');
      LatticePropagator ps = env.get_point(x, 's');
      // for(const auto &y: env.xgs_l)  {
        // peekSite(pl_site[x][y], pl, Coordinate(y));
        // peekSite(ps_site[x][y], ps, Coordinate(y));
      // }
      autoView(pl_v, pl, CpuRead);
      autoView(ps_v, ps, CpuRead);
      for(const auto &y: env.xgs_l)  {
        Coordinate ycoor = Coordinate(y);
        peekLocalSite(pl_site[x][y], pl_v, ycoor); // We used "assert to ensure that only one process is running
        peekLocalSite(ps_site[x][y], ps_v, ycoor);
      }
      num_pt_src++;
    }
    std::cout << GridLogMessage << "finished peeking" << std::endl;
    

    vector<double> avg_corr(T, 0.);  // correlation of each separation
    vector<int> counts(T, 0); // counts of pairs for each separation
    for(const auto &x: env.xgs_l) { // x is position of eta source
      for(const auto &y: env.xgs_l)  {

        LatticeComplexSite amp_connected = - 1./3. * trace(pl_site[x][y] * adj(pl_site[x][y])) - 2./3. * trace(ps_site[x][y] * adj(ps_site[x][y]));
        LatticeComplexSite amp_disconnected = 2./3. * ( trace(g5 * pl_site[x][x]) * trace(g5 * pl_site[y][y]) + trace(g5 * ps_site[x][x]) * trace(g5 * ps_site[y][y]) - trace(g5 * pl_site[x][x]) * trace(g5 * ps_site[y][y]) - trace(g5 * ps_site[x][x]) * trace(g5 * pl_site[y][y])  );
        double corr = ( amp_connected()()() + amp_disconnected()()() ).real();

        int tx = x[3], ty=y[3];
        int t_sep = (tx - ty + T) % T;
        avg_corr[t_sep] += corr;
        counts[t_sep] += 1;
      }
    }

    
    for(int t_sep=0; t_sep<T; ++t_sep) { // average over all point_sources
      avg_corr[t_sep] /= double(counts[t_sep]);
    }

    std::cout << GridLogMessage << "traj [" << traj << "]: " << avg_corr << std::endl;
  }

  std::cout << "Finished!" << std::endl;
  Grid_finalize();
  return 0;
}
