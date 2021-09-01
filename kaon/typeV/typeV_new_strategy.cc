#include "kaon/kaon.h"
#include "typeV_diagrams.h"
#include "../table3d_utils.h"

using namespace std;
using namespace Grid;

inline double calc_3d_vec_norm(int x, int y, int z) {
  return std::sqrt(double(x)*x + y*y + z*z);
}

// generate table4d[tK][xt][|vec{x}|][vt]; initialize it to 0
vector<vector<vector<vector<Complex>>>> initialize_typeV_table4d(int T, int num_R) { 
    vector<vector<vector<vector<Complex>>>> table4d(T);  
    for(int tK=0; tK<T; ++tK) {
      table4d[tK].resize(T);
      for(int xt=0; xt<T; ++xt) {
        table4d[tK][xt].resize(num_R);
        for(int R=0; R<num_R; ++R) {
          table4d[tK][xt][R].resize(T, 0.);
      }
    }
    return table4d;
}

vector<vector<LatticeComplexSite>> sum_typeV(const LatticeComplex &lat) {

  const int T = env.grid->_fdimensions[3];
  int X = env.grid->_fdimensions[0], Y = env.grid->_fdimensions[1], Z = env.grid->_fdimensions[2];
  int num_R = int(calc_3d_vec_norm(X-1, Y-1, Z-1));

  vector<vector<LatticeComplexSite>> rst(T); // rst[xt][R]
  for(int xt=0; xt<T; ++xt) rst[xt].resize(num_R, 0.); // initialize with 0.

  autoView(lat_v, lat, CpuRead);
  // thread_for(ss, lat.Grid()->lSites(), {
  for(int ss=0; ss<lat.Grid()->lSites(); ++ss) { // FIXME: cannot use thread for because we add to rst[t][int(R)]
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

		LatticeComplexSite m;
		peekLocalSite(m, lat_v, lcoor);

    const int x=gcoor[0], y=gcoor[1], z=gcoor[2], t=gcoor[3];
    double R = calc_3d_vec_norm(x, y, z);

    rst[t][int(R)] += m;
  // });
  }
  return rst;
}

vector<vector<Complex>> typeV_table4d_to_table2d(const vector<vector<vector<vector<<Complex>>>> &table4d, 
                                              int tsep_min=6, int tsep_max=14) { // table4d[tK][xt][R][vt] -> table2d[R][vt] with xt = 0
  int T = table4d.size();
  int num_R = table4d[0][0].size();

  // 1. average over xt
  vector<vector<vector<Complex>>> table3d(T); // table3d[tK][R][vt]
  for(auto &x: table3d) {
    x.resize(num_R);
    for(auto &y: x) {
      y.resize(T, 0.);
    }
  }
  for(int tK=0; tK<T; ++tK) {
    for(int R=0; R<num_R; ++R) { 
      for(int vt=0; vt<T; ++vt) {
        for(int xt=0; xt<T; ++xt) {
          table3d[tK][R][vt] += table4d[(tK+xt)%T][xt][R][(vt+xt)%T];
        }
        table3d[tK][R][vt] /= double(T);
      }
    }
  }

  // 2. average over tK
  vector<vector<Complex>> table2d(num_R); // table2d[R][vt]
  for(auto &x: table2d) x.resize(T, 0.);

  for(int R=0; R<num_R; ++R) { 
    for(int vt=0; vt<T; ++vt) {
      for(int tK=-max_tsep; tK<=-min_tsep; ++tK) { // average tK in [-max_TK, min_TK]
        table2d[R][vt] += table3d[(tK+T)%T][R][vt];
      }
      table2d[R][vt] /= double(max_tsep - min_tsep + 1);
    }
  }

  return table2d;
}



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
  map<string, decltype(&typeV)> amplitude_func {{"typeV", typeV}};

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
  int X = env.grid->_fdimensions[0], Y = env.grid->_fdimensions[1], Z = env.grid->_fdimensions[2];
  int num_R = int(calc_3d_vec_norm(X-1, Y-1, Z-1));

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    vector<vector<vector<vector<Complex>>>> table4d_D1Q1=initialize_typeV_table4d(T, num_R); // table4d[tK][xt][R][vt]
    vector<vector<vector<vector<Complex>>>> table4d_D1Q2=initialize_typeV_table4d(T, num_R);
    vector<vector<vector<vector<Complex>>>> table4d_D2Q1=initialize_typeV_table4d(T, num_R); 
    vector<vector<vector<vector<Complex>>>> table4d_D2Q2=initialize_typeV_table4d(T, num_R);
    vector<vector<vector<vector<Complex>>>> table4d_sBar_d_D1=initialize_typeV_table4d(T, num_R);
    vector<vector<vector<vector<Complex>>>> table4d_sBar_d_D2=initialize_typeV_table4d(T, num_R);

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

      LatticePropagator pl = env.get_point(v, 'l'); // pl = L(x, v) 
      LatticePropagator ps = env.get_point(v, 's'); // ps = H(x, v)

      for(int tK=0; tK<T; ++tK) {
        LatticeComplex rst_D1Q1(env.grid), rst_D1Q2(env.grid), rst_D2Q1(env.grid), rst_D2Q2(env.grid), sBar_d_D1(env.grid), sBar_d_D2(env.grid);
        amplitude_func[diagram](v, tK, rst_D1Q1, rst_D1Q2, rst_D2Q1, rst_D2Q2, sBar_d_D1, sBar_d_D2, env, wl, ws, pl, ps, Lxx, max_uv_sep);

        // vector<LatticeComplexSite> rst_D1Q1_xt, rst_D1Q2_xt, rst_D2Q1_xt, rst_D2Q2_xt, sBar_d_D1_xt, sBar_d_D2_xt;
        // sliceSum(rst_D1Q1, rst_D1Q1_xt, Tdir);
        // sliceSum(rst_D1Q2, rst_D1Q2_xt, Tdir);
        // sliceSum(rst_D2Q1, rst_D2Q1_xt, Tdir);
        // sliceSum(rst_D2Q2, rst_D2Q2_xt, Tdir);
        // sliceSum(sBar_d_D1, sBar_d_D1_xt, Tdir);
        // sliceSum(sBar_d_D2, sBar_d_D2_xt, Tdir);
        vector<vector<LatticeComplexSite>> rst_D1Q1_xt_R, rst_D1Q2_xt_R, rst_D2Q1_xt_R, rst_D2Q2_xt_R, sBar_d_D1_xt_R, sBar_d_D2_xt_R;
        rst_D1Q1_xt_R = sum_typeV(rst_D1Q1);
        rst_D1Q2_xt_R = sum_typeV(rst_D1Q2);
        rst_D2Q1_xt_R = sum_typeV(rst_D2Q1);
        rst_D2Q2_xt_R = sum_typeV(rst_D2Q2);
        sBar_d_D1_xt_R = sum_typeV(sBar_d_D1);
        sBar_d_D2_xt_R = sum_typeV(sBar_d_D2);
        for(int xt=0; xt<T; ++xt) {
          for(int R=0; R<num_R; ++R) {
            table4d_D1Q1[tK][xt][R][vt] += rst_D1Q1_xt[xt][R]()()() / double(pt_counts_slices[vt]); // divide by number of point sources on each time slice
            table4d_D1Q2[tK][xt][R][vt] += rst_D1Q2_xt[xt][R]()()() / double(pt_counts_slices[vt]);
            table4d_D2Q1[tK][xt][R][vt] += rst_D2Q1_xt[xt][R]()()() / double(pt_counts_slices[vt]); 
            table4d_D2Q2[tK][xt][R][vt] += rst_D2Q2_xt[xt][R]()()() / double(pt_counts_slices[vt]);
            table4d_sBar_d_D1[tK][xt][R][vt] += sBar_d_D1_xt[xt][R]()()() / double(pt_counts_slices[vt]);
            table4d_sBar_d_D2[tK][xt][R][vt] += sBar_d_D2_xt[xt][R]()()() / double(pt_counts_slices[vt]);
          }
        }
      }

    } // end of point source loop

    std::cout << "Total number of point sources: " << num_pt_src << std::endl;
    std::cout << "Number of point sources on each time slice: " << pt_counts_slices << std::endl;

    // std::cout << "traj [" << traj << "] amplitude Q1 by time slice: " << table3d_Q1 << std::endl;
    // std::cout << "traj [" << traj << "] amplitude Q2 by time slice: " << table3d_Q2 << std::endl;

    vector<vector<Complex>> table2d_D1Q1 = typeV_table4d_to_table2d(table4d_D1Q1);  // table2d[R][vt], with xt=0
    vector<vector<Complex>> table2d_D1Q2 = typeV_table4d_to_table2d(table4d_D1Q2);
    vector<vector<Complex>> table2d_D2Q1 = typeV_table4d_to_table2d(table4d_D2Q1);  // table2d[R][vt], with xt=0
    vector<vector<Complex>> table2d_D2Q2 = typeV_table4d_to_table2d(table4d_D2Q2);
    vector<vector<Complex>> table2d_sBar_d_D1 = typeV_table4d_to_table2d(table4d_sBar_d_D1);
    vector<vector<Complex>> table2d_sBar_d_D2 = typeV_table4d_to_table2d(table4d_sBar_d_D2);
    std::cout << "traj [" << traj << "] table2d_D1Q1[R][vt]: " << table2d_D1Q1 << std::endl;
    std::cout << "traj [" << traj << "] table2d_D1Q2[R][vt]: " << table2d_D1Q2 << std::endl;
    std::cout << "traj [" << traj << "] table2d_D2Q1[R][vt]: " << table2d_D2Q1 << std::endl;
    std::cout << "traj [" << traj << "] table2d_D2Q2[R][vt]: " << table2d_D2Q2 << std::endl;
    std::cout << "traj [" << traj << "] table2d sBar_d_D1[R][vt]: " << table2d_sBar_d_D1 << std::endl;
    std::cout << "traj [" << traj << "] table2d sBar_d_D2[R][vt]: " << table2d_sBar_d_D2 << std::endl;

  } // end of traj loop

  std::cout << GridLogMessage << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}
