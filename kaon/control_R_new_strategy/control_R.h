#pragma once

namespace Grid {


using namespace std;

inline double calc_3d_vec_norm(int x, int y, int z) {
  return std::sqrt(double(x)*x + y*y + z*z);
}

inline int get_relative_dist(int x, int x0, int X) { // return a value in (-X/2, X/2]
  int dist = x - x0;
  if(dist > X/2) dist -= X;
  return dist;
}

vector<vector<Complex>> sum_T_R(const LatticeComplex &lat, const std::vector<int> &v) { // v is a 3-element vector, center for calculating spatial distance

  // return a 2D vector rst[T][R]

  const int T = lat.Grid()->_fdimensions[3];
  int X = lat.Grid()->_fdimensions[0], Y = lat.Grid()->_fdimensions[1], Z = lat.Grid()->_fdimensions[2];
  // int num_R = 1 + int(calc_3d_vec_norm(X-1, Y-1, Z-1));
  int num_R = 1 + int(calc_3d_vec_norm(X/2, Y/2, Z/2));

  // std::cout << "X: " << X << " T: " << T << std::endl;
  // std::cout << "num_R: " << num_R << std::endl;
  vector<vector<Complex>> rst(T); // rst[xt][R]
  for(int xt=0; xt<T; ++xt) rst[xt].resize(num_R, 0.); // initialize with 0.

  autoView(lat_v, lat, CpuRead);
  // thread_for(ss, lat.Grid()->lSites(), {
  for(int ss=0; ss<lat.Grid()->lSites(); ++ss) { // FIXME: cannot use thread for because we add to rst[t][int(R)]
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

		LatticeComplexSite m;
		peekLocalSite(m, lat_v, lcoor);

    // int x = (gcoor[0] - v[0])
    int dist_x = get_relative_dist(gcoor[0], v[0], X);
    int dist_y = get_relative_dist(gcoor[1], v[1], Y);
    int dist_z = get_relative_dist(gcoor[2], v[2], Z);
    double R = calc_3d_vec_norm(dist_x, dist_y, dist_z);

    int t = gcoor[3];
    // std::cout << x << " " << y << " " << z << "t" << t << " int(R) " << int(R) << std::endl;
    rst[t][int(R)] += m()()();
  // });
  }

  for(int t=0; t<T; ++t) { // sum over all nodes
    double *ptr = (double *) rst[t].data();  
    int count = rst[t].size() * 2;
    lat.Grid()->GlobalSumVector(ptr, count); // Sum over Nodes
  }
  return rst;
}

// generate table4d[tK][xt][|vec{x}|][vt]; initialize it to 0
vector<vector<vector<vector<Complex>>>> control_R_initialize_table4d(const Coordinate &fdims) { 

    const int T = fdims[3];
    int X = fdims[0], Y = fdims[1], Z = fdims[2];
    // int num_R = 1 + int(calc_3d_vec_norm(X-1, Y-1, Z-1));
    int num_R = 1 + int(calc_3d_vec_norm(X/2, Y/2, Z/2));

    vector<vector<vector<vector<Complex>>>> table4d(T);  
    for(int tK=0; tK<T; ++tK) {
      table4d[tK].resize(T);
      for(int xt=0; xt<T; ++xt) {
        table4d[tK][xt].resize(num_R);
        for(int R=0; R<num_R; ++R) {
          table4d[tK][xt][R].resize(T, 0.);
        }
      }
    }
    return table4d;
}

// input table4d[T][xt][R][vt], average T and (xt+vt) to get table2d[R][vt] with xt = 0
// table4d[tK][xt][R][vt] -> table2d[R][vt] with xt = 0
vector<vector<Complex>> control_R_table4d_to_table2d(const vector<vector<vector<vector<Complex>>>> &table4d, 
                                              int min_tsep=6, int max_tsep=14) { 
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
  // std::cout << "after averging over xt" << std::endl;

  // 2. average over tK
  vector<vector<Complex>> table2d(num_R); // table2d[R][vt]
  for(auto &x: table2d) x.resize(T, 0.);

  for(int R=0; R<num_R; ++R) { 
    for(int vt=0; vt<T; ++vt) {
      for(int tK=-max_tsep; tK<=-min_tsep; ++tK) { // average tK in [-max_TK, -min_TK]
        // std::cout << R << " " << vt << " " << tK << " " << (tK+T)%T << std::endl;
        // std::cout << table3d.size() << " " << table3d[0].size() << std::endl;
        // std::cout << table3d[(tK+T)%T][R][vt] << std::endl;
        // std::cout << table2d[R][vt] << std::endl;
        table2d[R][vt] += table3d[(tK+T)%T][R][vt];
      }
      table2d[R][vt] /= double(max_tsep - min_tsep + 1);
    }
  }

  return table2d;
}

}
