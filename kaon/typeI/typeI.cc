#include "typeI_utils.h"
#include "typeI_diagrams.h"

using namespace std;
using namespace Grid;

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);
  // zyd_init_Grid_Qlattice(argc, argv);
 
#ifdef CUTH_FREE_FIELD
  std::vector<int> tseps = {1,2};
  int max_uv_sep = 3;
  Env env("FreeField_8nt8");
#else
  std::vector<int> tseps = {6,8,10,12,14}; // FIXME
  int max_uv_sep = 16;
  Env env("24ID");
#endif

  map<string, decltype(&typeI_D1a)> amplitude_func {{"typeI_D1a", typeI_D1a},  {"typeI_D1b", typeI_D1b}, 
                                                    {"typeI_D2a", typeI_D2a},  {"typeI_D2b", typeI_D2b}};

  // env.N_pt_src = 1;  // keep only one point
  env.N_pt_src = -1;  // Use all points

  /////////////////////////////////////////////////

  const int T = env.grid->_fdimensions[3];

  std::cout << "tseps: " << tseps << std::endl;
  std::cout << "max_uv_sep: " << max_uv_sep << std::endl;

  string diagram;
  if( GridCmdOptionExists(argv, argv+argc, "--diagram") ) {
    diagram = GridCmdOptionPayload(argv, argv+argc, "--diagram");
    std::cout << "Calculating diagram " << diagram << std::endl;
    assert(amplitude_func.find(diagram) != amplitude_func.end());
  }
  else {
    std::cout << "--diagram not specified; exiting" << std::endl;
    exit(0);
  }

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

  // int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    vector<vector<Complex>> table2d_Q1(tseps.size()), table2d_Q2(tseps.size());
    for(int i=0; i<tseps.size(); ++i) {
      table2d_Q1[i].resize(T, 0.);
      table2d_Q2[i].resize(T, 0.);
    }

    vector<LatticePropagator> wl = env.get_wall('l');
    vector<LatticePropagator> ws = env.get_wall('s');

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_s.resize(env.N_pt_src);
    for(const auto &x: env.xgs_s) {
      std::cout << "# Point source: " << num_pt_src << std::endl;
      std::cout << "Point source: " << x << std::endl;
      // print_memory();

      ++num_pt_src;

      LatticePropagator pl = env.get_point(x, 'l'); 
      LatticePropagator ps = env.get_point(x, 's'); 

      for(int tsep_idx=0; tsep_idx<tseps.size(); ++tsep_idx) {
        int tsep = tseps[tsep_idx];
        int tK = (x[3] - tsep + T) % T;

        LatticeComplex rst_Q1(env.grid), rst_Q2(env.grid);
        amplitude_func[diagram](x, tK, rst_Q1, rst_Q2, env, wl, ws, pl, ps, max_uv_sep);

        vector<LatticeComplexSite> rst_Q1_vt, rst_Q2_vt; // Sum over each time slice of v
        sliceSum(rst_Q1, rst_Q1_vt, Tdir);
        sliceSum(rst_Q2, rst_Q2_vt, Tdir);

        for(int vt=0; vt<T; ++vt) {
          table2d_Q1[tsep_idx][vt] += rst_Q1_vt[(vt+x[3])%T]()()(); // For table2d, xt := 0
          table2d_Q2[tsep_idx][vt] += rst_Q2_vt[(vt+x[3])%T]()()();
        }
      }

    } // end of point source loop

    std::cout << GridLogMessage << "Number of point sources: " << num_pt_src << std::endl;
    for(int tsep_idx=0; tsep_idx<tseps.size(); ++tsep_idx) {
      for(int i=0; i<T; ++i) {
        table2d_Q1[tsep_idx][i] /= double(num_pt_src);
        table2d_Q2[tsep_idx][i] /= double(num_pt_src);
      }
    }

    std::cout << "traj [" << traj << "] table2d_Q1: " << table2d_Q1 << std::endl;
    std::cout << "traj [" << traj << "] table2d_Q2: " << table2d_Q2 << std::endl;

  } // end of traj loop

  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}
