#include "kaon/kaon.h"
#include "typeIV_a_diagrams.h"
#include "../table3d_utils.h"

using namespace std;
using namespace Grid;


int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

#ifdef CUTH_FREE_FIELD
  int max_uv_sep = 3;
  Env env("FreeField_8nt8");
#else
  int max_uv_sep = 16;
  Env env("24ID");
#endif

  env.N_pt_src = -1;  
  map<string, decltype(&typeIV_D1a)> amplitude_func {{"typeIV_D1a", typeIV_D1a},  {"typeIV_D2a", typeIV_D2a}};

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

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    vector<vector<vector<Complex>>> table3d_Q1=initialize_table3d(T); // table3d[tK][xt][vt]
    vector<vector<vector<Complex>>> table3d_Q2=initialize_table3d(T);
    vector<vector<vector<Complex>>> table3d_sBar_d=initialize_table3d(T);

    std::vector<LatticePropagator> wl = env.get_wall('l');
    std::vector<LatticePropagator> ws = env.get_wall('s');
    LatticePropagator Lxx = env.get_Lxx();

    if(env.N_pt_src != -1) env.xgs_l.resize(env.N_pt_src);

    std::vector<int> pt_counts_slices(T, 0);
    for(const auto &pt: env.xgs_l) ++pt_counts_slices[pt[3]]; // Number of point sources on each time slice
    for(int i=0; i<T; ++i) assert(pt_counts_slices[i] > 0); // Must have at least one point source on each time slice

    int num_pt_src = 0;
    for(const auto &v: env.xgs_l) {
      std::cout << "# Point source: " << num_pt_src << std::endl;
      std::cout << "Point source: " << v << std::endl;
      ++num_pt_src;

      int vt = v[3];

      LatticePropagator point_prop(env.grid);
      if(diagram=="typeIV_D1a") point_prop = env.get_point(v, 'l'); // pl = L(x, v) // for typeIV_D1a 
      else point_prop = env.get_point(v, 's'); // ps = H(x, v) // for typeIV_D2a

      for(int tK=0; tK<T; ++tK) {
        LatticeComplex rst_Q1(env.grid), rst_Q2(env.grid), rst_sBar_d(env.grid);
        amplitude_func[diagram](v, tK, rst_Q1, rst_Q2, rst_sBar_d, env, wl, ws, point_prop, Lxx, max_uv_sep);

        vector<LatticeComplexSite> rst_Q1_xt, rst_Q2_xt, rst_sBar_d_xt;
        sliceSum(rst_Q1, rst_Q1_xt, Tdir);
        sliceSum(rst_Q2, rst_Q2_xt, Tdir);
        sliceSum(rst_sBar_d, rst_sBar_d_xt, Tdir);
        for(int xt=0; xt<T; ++xt) {
          table3d_Q1[tK][xt][vt] += rst_Q1_xt[xt]()()() / double(pt_counts_slices[vt]); // divide by number of point sources on each time slice
          table3d_Q2[tK][xt][vt] += rst_Q2_xt[xt]()()() / double(pt_counts_slices[vt]);
          table3d_sBar_d[tK][xt][vt] += rst_sBar_d_xt[xt]()()() / double(pt_counts_slices[vt]);
        }
      }

    } // end of point source loop

    std::cout << "Total number of point sources: " << num_pt_src << std::endl;
    std::cout << "Number of point sources on each time slice: " << pt_counts_slices << std::endl;

    // std::cout << "traj [" << traj << "] amplitude Q1 by time slice: " << table3d_Q1 << std::endl;
    // std::cout << "traj [" << traj << "] amplitude Q2 by time slice: " << table3d_Q2 << std::endl;

    vector<vector<Complex>> table2d_Q1 = table3d_to_table2d(table3d_Q1);  // table2d[tK][vt], with xt=0
    vector<vector<Complex>> table2d_Q2 = table3d_to_table2d(table3d_Q2);
    vector<vector<Complex>> table2d_sBar_d = table3d_to_table2d(table3d_sBar_d);
    std::cout << "traj [" << traj << "] amplitude Q1: " << table2d_Q1 << std::endl;
    std::cout << "traj [" << traj << "] amplitude Q2: " << table2d_Q2 << std::endl;
    std::cout << "traj [" << traj << "] amplitude sBar_d T2 diagram a: " << table2d_sBar_d << std::endl;

  } // end of traj loop

  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}
